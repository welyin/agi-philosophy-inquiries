#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""手动做 git gc：删除不可到达的 objects"""

import os
import glob
from dulwich.repo import Repo

def main():
    repo_path = r'D:\workspace\AGI的哲学思考'
    repo = Repo(repo_path)
    
    print("Collecting all reachable objects...")
    
    # 从所有引用开始，遍历所有可到达的 objects
    reachable = set()
    
    def add_object(sha):
        if sha in reachable:
            return
        reachable.add(sha)
        obj = repo.object_store[sha]
        
        if isinstance(obj, type(repo[repo.refs[b'refs/heads/main']])):  # Commit
            add_object(obj.tree)
            for parent in obj.parents:
                add_object(parent)
        elif hasattr(obj, 'items'):  # Tree
            for entry in obj.items():
                add_object(entry.sha)
        # Blob 不需要递归
    
    # 遍历所有引用
    for ref_name, ref_sha in repo.get_refs().items():
        ref_name_str = ref_name.decode()
        if ref_name_str.startswith('refs/codex/'):
            continue  # 跳过 codex 的临时引用
        print(f"  Following ref: {ref_name_str}")
        add_object(ref_sha)
    
    print(f"\nReachable objects: {len(reachable)}")
    
    # 收集所有 objects
    objects_dir = os.path.join(repo_path, '.git', 'objects')
    all_objects = set()
    
    for dirpath, dirnames, filenames in os.walk(objects_dir):
        # 跳过 pack 目录（先不处理 packfile）
        if 'pack' in dirpath:
            continue
        for fname in filenames:
            if len(fname) == 38:  # sha1 是 40 字符，前 2 字符是目录名
                sha = os.path.basename(dirpath) + fname
                all_objects.add(sha.encode())
    
    print(f"Total objects: {len(all_objects)}")
    
    # 计算不可到达的 objects
    unreachable = all_objects - reachable
    print(f"Unreachable objects: {len(unreachable)}")
    
    # 删除不可到达的 objects
    deleted = 0
    for sha in unreachable:
        sha_str = sha.decode()
        path = os.path.join(objects_dir, sha_str[:2], sha_str[2:])
        if os.path.exists(path):
            try:
                os.remove(path)
                deleted += 1
            except Exception as e:
                pass
    
    print(f"Deleted {deleted} unreachable objects")
    
    # 清理空目录
    for dirpath, dirnames, filenames in os.walk(objects_dir, topdown=False):
        if dirpath == objects_dir:
            continue
        try:
            if not os.listdir(dirpath):
                os.rmdir(dirpath)
        except:
            pass
    
    print("Done.")

if __name__ == "__main__":
    main()
