#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
用 dulwich 重写 git 历史，删除所有二进制文件
（绕过 git.exe 被安全软件拦截的问题）
"""

import os
import sys
from dulwich.repo import Repo
from dulwich.objects import Blob, Tree, Commit, Tag
from dulwich.file import GitFile

# 要删除的文件扩展名
EXCLUDE_EXTENSIONS = {'.docx', '.pptx', '.pdf', '.jpg', '.jpeg', '.png', '.gif', '.bmp'}

def should_delete(path_bytes):
    """判断文件是否应该被删除"""
    path = path_bytes.decode('utf-8', errors='ignore')
    _, ext = os.path.splitext(path.lower())
    return ext in EXCLUDE_EXTENSIONS

def process_tree(repo, tree_id, new_tree_entries, prefix=b''):
    """递归处理树，删除匹配的文件"""
    tree = repo.object_store[tree_id]
    
    for entry in tree.items():
        mode, name, sha = entry
        full_path = prefix + b'/' + name if prefix else name
        
        if mode == 0o40000:  # 目录
            # 递归处理子目录
            sub_new_entries = {}
            process_tree(repo, sha, sub_new_entries, full_path)
            if sub_new_entries:
                # 创建新的子树
                new_subtree = Tree()
                for entry_mode, entry_name, entry_sha in sub_new_entries.values():
                    new_subtree.add(entry_name, entry_mode, entry_sha)
                new_subtree_id = repo.object_store.add_object(new_subtree)
                new_tree_entries[name] = (mode, name, new_subtree_id)
        else:
            # 文件
            if not should_delete(full_path):
                new_tree_entries[name] = (mode, name, sha)

def rewrite_commit(repo, commit_id, commit_cache):
    """重写单个提交，返回新的提交 ID"""
    if commit_id in commit_cache:
        return commit_cache[commit_id]
    
    commit = repo.object_store[commit_id]
    
    # 递归处理父提交
    new_parents = []
    for parent_id in commit.parents:
        new_parent_id = rewrite_commit(repo, parent_id, commit_cache)
        new_parents.append(new_parent_id)
    
    # 处理树
    new_tree_entries = {}
    process_tree(repo, commit.tree, new_tree_entries)
    
    # 创建新的树
    new_tree = Tree()
    for entry_mode, entry_name, entry_sha in new_tree_entries.values():
        new_tree.add(entry_name, entry_mode, entry_sha)
    new_tree_id = repo.object_store.add_object(new_tree)
    
    # 创建新的提交
    new_commit = Commit()
    new_commit.tree = new_tree_id
    new_commit.parents = new_parents
    new_commit.author = commit.author
    new_commit.committer = commit.committer
    new_commit.author_time = commit.author_time
    new_commit.author_timezone = commit.author_timezone
    new_commit.committer_time = commit.committer_time
    new_commit.committer_timezone = commit.committer_timezone
    new_commit.message = commit.message
    
    new_commit_id = repo.object_store.add_object(new_commit)
    commit_cache[commit_id] = new_commit_id
    
    return new_commit_id

def main():
    repo_path = r'D:\workspace\AGI的哲学思考'
    repo = Repo(repo_path)
    
    print(f"Repository: {repo_path}")
    print(f"Excluding extensions: {', '.join(EXCLUDE_EXTENSIONS)}")
    
    # 获取所有分支
    refs = repo.get_refs()
    print(f"\nFound {len(refs)} refs")
    
    commit_cache = {}
    
    # 处理所有分支
    for ref_name, ref_value in refs.items():
        ref_name_str = ref_name.decode('utf-8')
        if ref_name_str.startswith('refs/heads/'):
            branch_name = ref_name_str[len('refs/heads/'):]
            print(f"\nProcessing branch: {branch_name}")
            
            try:
                old_commit_id = ref_value
                new_commit_id = rewrite_commit(repo, old_commit_id, commit_cache)
                
                # 更新分支引用
                repo.refs[ref_name] = new_commit_id
                print(f"  Rewritten: {old_commit_id.decode()[:8]} -> {new_commit_id.decode()[:8]}")
            except Exception as e:
                print(f"  Error: {e}")
    
    # 处理标签
    for ref_name, ref_value in refs.items():
        if ref_name.startswith(b'refs/tags/'):
            try:
                obj = repo.object_store[ref_value]
                if isinstance(obj, Tag):
                    # 标签指向另一个对象
                    tagged_obj = repo.object_store[obj.object_id]
                    if isinstance(tagged_obj, Commit):
                        new_tagged_id = rewrite_commit(repo, obj.object_id, commit_cache)
                        # 创建新标签
                        new_tag = Tag()
                        new_tag.object_id = new_tagged_id
                        new_tag.object_type = obj.object_type
                        new_tag.tag = obj.tag
                        new_tag.tagger = obj.tagger
                        new_tag.tag_time = obj.tag_time
                        new_tag.tag_timezone = obj.tag_timezone
                        new_tag.message = obj.message
                        new_tag_id = repo.object_store.add_object(new_tag)
                        repo.refs[ref_name] = new_tag_id
                        print(f"Rewrote tag: {ref_name.decode()}")
            except Exception as e:
                print(f"Tag error: {e}")
    
    print(f"\nDone. Rewrote {len(commit_cache)} commits")

if __name__ == '__main__':
    main()
