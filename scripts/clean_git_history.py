#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""从 git 历史中删除所有二进制文件（docx, pptx, pdf, jpg, png 等）"""

import os
import re
from dulwich.repo import Repo
from dulwich.objects import Commit, Tree, Blob

# 要删除的文件扩展名
PATTERNS = ['.docx', '.pptx', '.pdf', '.jpg', '.jpeg', '.png', '.gif', '.bmp']

def should_delete(path):
    """判断文件是否应该被删除"""
    lower = path.lower()
    for pat in PATTERNS:
        if lower.endswith(pat):
            return True
    return False

def filter_tree(repo, tree_id, prefix=""):
    """递归过滤树，返回新的 tree_id"""
    tree = repo.object_store[tree_id]
    
    new_entries = []
    deleted = []
    
    for entry in tree.items():
        name = entry.path.decode() if isinstance(entry.path, bytes) else entry.path
        full_path = prefix + name
        
        mode = entry.mode
        sha = entry.sha
        
        if mode == 0o40000:  # 目录
            # 递归处理子目录
            new_subtree = filter_tree(repo, sha, full_path + "/")
            if new_subtree is not None:
                new_entries.append((mode, entry.path, new_subtree))
        else:  # 文件
            if should_delete(full_path):
                deleted.append(full_path)
            else:
                new_entries.append((mode, entry.path, sha))
    
    if not new_entries:
        # 空目录返回 None
        return None
    
    # 创建新的树
    new_tree = Tree()
    for mode, path, sha in new_entries:
        new_tree.add(path, mode, sha)
    
    repo.object_store.add_object(new_tree)
    return new_tree.id

def filter_commit(repo, commit_id, new_parents_cache):
    """递归过滤提交，返回新的 commit_id"""
    commit = repo.object_store[commit_id]
    
    # 先处理父提交
    new_parents = []
    for parent_id in commit.parents:
        if parent_id in new_parents_cache:
            new_parents.append(new_parents_cache[parent_id])
        else:
            new_parent = filter_commit(repo, parent_id, new_parents_cache)
            new_parents.append(new_parent)
            new_parents_cache[parent_id] = new_parent
    
    # 过滤树
    new_tree_id = filter_tree(repo, commit.tree)
    
    # 创建新的提交
    new_commit = Commit()
    new_commit.tree = new_tree_id
    new_commit.parents = new_parents
    new_commit.author = commit.author
    new_commit.committer = commit.committer
    new_commit.commit_time = commit.commit_time
    new_commit.commit_timezone = commit.commit_timezone
    new_commit.author_time = commit.author_time
    new_commit.author_timezone = commit.author_timezone
    new_commit.message = commit.message
    
    repo.object_store.add_object(new_commit)
    
    new_parents_cache[commit_id] = new_commit.id
    return new_commit.id

def main():
    repo_path = r'D:\workspace\AGI的哲学思考'
    repo = Repo(repo_path)
    
    print(f"Repo: {repo_path}")
    print(f"Deleting files matching: {PATTERNS}")
    print()
    
    # 找到所有分支引用
    refs = repo.get_refs()
    print(f"Found {len(refs)} refs")
    
    # 缓存：旧 commit_id -> 新 commit_id
    commit_cache = {}
    
    # 处理每个分支
    for ref_name, ref_sha in refs.items():
        ref_name_str = ref_name.decode() if isinstance(ref_name, bytes) else ref_name
        print(f"\nProcessing ref: {ref_name_str}")
        
        # 只处理分支和标签，不处理临时的
        if ref_name_str.startswith('refs/codex/'):
            print(f"  Skipping: {ref_name_str}")
            continue
        
        try:
            obj = repo.object_store[ref_sha]
            if isinstance(obj, Commit):
                new_sha = filter_commit(repo, ref_sha, commit_cache)
                print(f"  Old: {ref_sha.decode()[:8]}")
                print(f"  New: {new_sha.decode()[:8]}")
                
                # 更新引用
                repo.refs[ref_name] = new_sha
            else:
                print(f"  Not a commit: {type(obj)}")
        except Exception as e:
            print(f"  Error: {e}")
    
    # 清理旧的 objects（可选，先不做）
    print(f"\nDone. Modified {len(commit_cache)} commits.")
    print("Now run: git reflog expire --expire=now --all && git gc --prune=now --aggressive")

if __name__ == "__main__":
    main()
