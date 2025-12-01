#!/usr/bin/env python3
"""
Script to recursively find files with non-ASCII characters in their names.

This is a belt and suspenders routine, in the sense that it just checks 
for non_ascii filenames, see rename_unicode_filenames
"""

import os
import sys
import argparse
from pathlib import Path


def is_ascii(text):
    """Check if a string contains only ASCII characters."""
    try:
        text.encode('ascii')
        return True
    except UnicodeEncodeError:
        return False


def find_non_ascii_files(directory, show_ascii=False):
    """
    Recursively find files with non-ASCII names in the given directory.
    
    Args:
        directory (str): Path to the directory to scan
        show_ascii (bool): If True, also show ASCII filenames for comparison
    
    Returns:
        tuple: (non_ascii_files, ascii_files) - lists of file paths
    """
    non_ascii_files = []
    ascii_files = []
    
    try:
        for root, dirs, files in os.walk(directory):
            # Check files
            for file in files:
                full_path = os.path.join(root, file)
                if is_ascii(file):
                    ascii_files.append(full_path)
                else:
                    non_ascii_files.append(full_path)
            
            # Check directory names too
            for dir_name in dirs:
                if not is_ascii(dir_name):
                    full_path = os.path.join(root, dir_name)
                    non_ascii_files.append(full_path + "/")  # Mark as directory
    
    except PermissionError as e:
        print(f"Permission denied: {e}", file=sys.stderr)
    except Exception as e:
        print(f"Error scanning directory: {e}", file=sys.stderr)
    
    return non_ascii_files, ascii_files


def main():
    parser = argparse.ArgumentParser(
        description="Find files and directories with non-ASCII characters in their names"
    )
    parser.add_argument(
        "directory", 
        nargs="?", 
        default=".", 
        help="Directory to scan (default: current directory)"
    )
    parser.add_argument(
        "--show-ascii", 
        action="store_true", 
        help="Also display files with ASCII names"
    )
    parser.add_argument(
        "--count-only", 
        action="store_true", 
        help="Only show counts, not file names"
    )
    
    args = parser.parse_args()
    
    # Validate directory
    if not os.path.isdir(args.directory):
        print(f"Error: '{args.directory}' is not a valid directory", file=sys.stderr)
        sys.exit(1)
    
    print(f"Scanning directory: {os.path.abspath(args.directory)}")
    print("=" * 50)
    
    non_ascii_files, ascii_files = find_non_ascii_files(args.directory, args.show_ascii)
    
    if args.count_only:
        print(f"Files/directories with non-ASCII names: {len(non_ascii_files)}")
        if args.show_ascii:
            print(f"Files/directories with ASCII names: {len(ascii_files)}")
    else:
        if non_ascii_files:
            print(f"Found {len(non_ascii_files)} files/directories with non-ASCII names:")
            print("-" * 40)
            for file_path in sorted(non_ascii_files):
                # Show the problematic characters more clearly
                filename = os.path.basename(file_path.rstrip('/'))
                try:
                    # Try to represent non-ASCII chars
                    print(f"{file_path}")
                    print(f"  → Filename bytes: {filename.encode('utf-8')}")
                except Exception as e:
                    print(f"{file_path} (encoding error: {e})")
                print()
        else:
            print("No files or directories with non-ASCII names found.")
        
        if args.show_ascii and ascii_files:
            print(f"\nFound {len(ascii_files)} files/directories with ASCII names:")
            print("-" * 40)
            for file_path in sorted(ascii_files):
                print(file_path)
    
    # Exit with non-zero code if non-ASCII files were found
    sys.exit(0 if not non_ascii_files else 1)


if __name__ == "__main__":
    main()
