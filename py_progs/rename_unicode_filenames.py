#!/usr/bin/env python3
"""
Script to recursively find files with non-ASCII characters in their names
and rename them to ASCII equivalents.
"""

import os
import sys
import argparse
import unicodedata
import re
from pathlib import Path
import shutil


# Custom mapping for filename-safe Unicode characters to ASCII equivalents
FILENAME_UNICODE_TO_ASCII_MAP = {
    # Quotation marks (often problematic in filenames)
    '\u2018': "'",      # LEFT SINGLE QUOTATION MARK
    '\u2019': "'",      # RIGHT SINGLE QUOTATION MARK
    '\u201A': "'",      # SINGLE LOW-9 QUOTATION MARK
    '\u201B': "'",      # SINGLE HIGH-REVERSED-9 QUOTATION MARK
    '\u201C': '"',      # LEFT DOUBLE QUOTATION MARK
    '\u201D': '"',      # RIGHT DOUBLE QUOTATION MARK
    '\u201E': '"',      # DOUBLE LOW-9 QUOTATION MARK
    '\u201F': '"',      # DOUBLE HIGH-REVERSED-9 QUOTATION MARK
    '\u2039': "'",      # SINGLE LEFT-POINTING ANGLE QUOTATION MARK
    '\u203A': "'",      # SINGLE RIGHT-POINTING ANGLE QUOTATION MARK
    '\u00AB': '"',      # LEFT-POINTING DOUBLE ANGLE QUOTATION MARK
    '\u00BB': '"',      # RIGHT-POINTING DOUBLE ANGLE QUOTATION MARK
    
    # Dashes and hyphens
    '\u2010': '-',      # HYPHEN
    '\u2011': '-',      # NON-BREAKING HYPHEN
    '\u2012': '-',      # FIGURE DASH
    '\u2013': '-',      # EN DASH
    '\u2014': '-',      # EM DASH
    '\u2015': '-',      # HORIZONTAL BAR
    '\u2212': '-',      # MINUS SIGN
    
    # Spaces (convert special spaces to regular spaces)
    '\u00A0': ' ',      # NON-BREAKING SPACE
    '\u2000': ' ',      # EN QUAD
    '\u2001': ' ',      # EM QUAD
    '\u2002': ' ',      # EN SPACE
    '\u2003': ' ',      # EM SPACE
    '\u2004': ' ',      # THREE-PER-EM SPACE
    '\u2005': ' ',      # FOUR-PER-EM SPACE
    '\u2006': ' ',      # SIX-PER-EM SPACE
    '\u2007': ' ',      # FIGURE SPACE
    '\u2008': ' ',      # PUNCTUATION SPACE
    '\u2009': ' ',      # THIN SPACE
    '\u200A': ' ',      # HAIR SPACE
    '\u202F': ' ',      # NARROW NO-BREAK SPACE
    '\u205F': ' ',      # MEDIUM MATHEMATICAL SPACE
    '\u3000': ' ',      # IDEOGRAPHIC SPACE
    
    # Dots and bullets (convert to safe characters)
    '\u2022': '_',      # BULLET
    '\u2023': '_',      # TRIANGULAR BULLET
    '\u2024': '.',      # ONE DOT LEADER
    '\u2025': '..',     # TWO DOT LEADER
    '\u2026': '...',    # HORIZONTAL ELLIPSIS
    '\u00B7': '_',      # MIDDLE DOT
    '\u2219': '_',      # BULLET OPERATOR
    
    # Mathematical symbols (filename-safe versions)
    '\u00B1': 'pm',     # PLUS-MINUS SIGN
    '\u00D7': 'x',      # MULTIPLICATION SIGN
    '\u00F7': 'div',    # DIVISION SIGN
    '\u2260': 'ne',     # NOT EQUAL TO
    '\u2264': 'le',     # LESS-THAN OR EQUAL TO
    '\u2265': 'ge',     # GREATER-THAN OR EQUAL TO
    '\u00B0': 'deg',    # DEGREE SIGN
    '\u2032': "'",      # PRIME
    '\u2033': "''",     # DOUBLE PRIME
    '\u221E': 'inf',    # INFINITY
    
    # Currency symbols
    '\u00A2': 'c',      # CENT SIGN
    '\u00A3': 'GBP',    # POUND SIGN
    '\u00A4': '$',      # GENERIC CURRENCY SYMBOL (though $ might be problematic)
    '\u00A5': 'JPY',    # YEN SIGN
    '\u20AC': 'EUR',    # EURO SIGN
    
    # Miscellaneous symbols
    '\u00A9': '(C)',    # COPYRIGHT SIGN
    '\u00AE': '(R)',    # REGISTERED SIGN
    '\u2122': 'TM',     # TRADE MARK SIGN
    '\u00A7': 'S',      # SECTION SIGN
    '\u00B6': 'P',      # PILCROW SIGN
    '\u2020': '+',      # DAGGER
    '\u2021': '++',     # DOUBLE DAGGER
}


def has_non_ascii_chars(text):
    """Check if text contains any non-ASCII characters."""
    return any(ord(char) > 127 for char in text)


def remove_accents(text):
    """Remove accents from characters using Unicode normalization."""
    # Normalize to NFD (decomposed form) and filter out combining characters
    normalized = unicodedata.normalize('NFD', text)
    ascii_text = ''.join(
        char for char in normalized 
        if unicodedata.category(char) != 'Mn'  # Mn = Nonspacing_Mark (accents)
    )
    return ascii_text


def clean_filename_chars(filename):
    """Remove or replace characters that are problematic in filenames."""
    # Characters that are often problematic in filenames across different OS
    problematic_chars = {
        '<': '(',
        '>': ')',
        ':': '-',
        '|': '-',
        '?': '',
        '*': '_',
        '/': '-',
        '\\': '-',
        '"': "'",
    }
    
    result = filename
    for char, replacement in problematic_chars.items():
        result = result.replace(char, replacement)
    
    return result


def convert_filename_to_ascii(filename, 
                             use_custom_map=True, 
                             remove_accents_flag=True,
                             fallback_char='_',
                             clean_problematic=True):
    """
    Convert a filename with Unicode characters to ASCII equivalent.
    
    Args:
        filename (str): Original filename
        use_custom_map (bool): Use custom Unicode-to-ASCII mapping
        remove_accents_flag (bool): Remove accents from characters
        fallback_char (str): Character to use when no conversion is possible
        clean_problematic (bool): Clean characters that are problematic in filenames
    
    Returns:
        str: ASCII-safe filename
    """
    result = filename
    
    # Step 1: Apply custom mapping
    if use_custom_map:
        for unicode_char, ascii_equiv in FILENAME_UNICODE_TO_ASCII_MAP.items():
            result = result.replace(unicode_char, ascii_equiv)
    
    # Step 2: Remove accents (decompose and remove combining marks)
    if remove_accents_flag:
        result = remove_accents(result)
    
    # Step 3: Handle remaining non-ASCII characters
    ascii_result = []
    for char in result:
        if ord(char) <= 127:
            # ASCII character, keep as-is
            ascii_result.append(char)
        else:
            # Try to find ASCII representation through normalization
            try:
                normalized = unicodedata.normalize('NFKD', char)
                ascii_version = ''.join(
                    c for c in normalized 
                    if ord(c) <= 127 or unicodedata.category(c) != 'Mn'
                )
                if ascii_version and all(ord(c) <= 127 for c in ascii_version):
                    ascii_result.append(ascii_version)
                else:
                    ascii_result.append(fallback_char)
            except:
                ascii_result.append(fallback_char)
    
    result = ''.join(ascii_result)
    
    # Step 4: Clean problematic characters
    if clean_problematic:
        result = clean_filename_chars(result)
    
    # Step 5: Clean up multiple spaces/underscores and trim
    result = re.sub(r'[ _]+', ' ', result)  # Multiple spaces/underscores to single space
    result = result.strip(' ._-')  # Remove leading/trailing problematic chars
    
    # Step 6: Ensure we don't have an empty filename
    if not result or result.isspace():
        result = 'renamed_file'
    
    return result


def generate_unique_filename(directory, base_name, extension):
    """
    Generate a unique filename in the given directory.
    
    Args:
        directory (str): Directory path
        base_name (str): Base filename (without extension)
        extension (str): File extension (including dot)
    
    Returns:
        str: Unique filename
    """
    counter = 1
    new_name = base_name + extension
    
    while os.path.exists(os.path.join(directory, new_name)):
        new_name = f"{base_name}_{counter}{extension}"
        counter += 1
    
    return new_name


def rename_file(old_path, new_filename, dry_run=False):
    """
    Rename a file to the new filename.
    
    Args:
        old_path (str): Current file path
        new_filename (str): New filename (not full path)
        dry_run (bool): If True, don't actually rename
    
    Returns:
        dict: Rename operation results
    """
    directory = os.path.dirname(old_path)
    old_filename = os.path.basename(old_path)
    
    # Handle extension
    old_name, old_ext = os.path.splitext(old_filename)
    new_name, new_ext = os.path.splitext(new_filename)
    
    # If new filename doesn't have extension, use original extension
    if not new_ext and old_ext:
        new_filename = new_name + old_ext
    
    # Ensure filename is unique
    final_filename = generate_unique_filename(directory, 
                                            os.path.splitext(new_filename)[0], 
                                            os.path.splitext(new_filename)[1])
    
    new_path = os.path.join(directory, final_filename)
    
    result = {
        'old_path': old_path,
        'new_path': new_path,
        'old_filename': old_filename,
        'new_filename': final_filename,
        'success': False,
        'error': None,
        'was_renamed': final_filename != old_filename
    }
    
    if not dry_run and result['was_renamed']:
        try:
            os.rename(old_path, new_path)
            result['success'] = True
        except Exception as e:
            result['error'] = str(e)
    else:
        result['success'] = True  # Dry run or no rename needed always "succeeds"
    
    return result


def find_and_rename_files(directory, 
                         rename_dirs=False, 
                         dry_run=False, 
                         **conversion_options):
    """
    Recursively find and rename files with non-ASCII characters in their names.
    
    Args:
        directory (str): Root directory to search
        rename_dirs (bool): Whether to also rename directories
        dry_run (bool): If True, don't actually rename files
        **conversion_options: Options for filename conversion
    
    Returns:
        dict: Results summary
    """
    results = {
        'total_files': 0,
        'total_dirs': 0,
        'files_with_non_ascii': 0,
        'dirs_with_non_ascii': 0,
        'files_renamed': 0,
        'dirs_renamed': 0,
        'files_processed': [],
        'dirs_processed': [],
        'errors': []
    }
    
    # We need to process directories from deepest to shallowest to avoid
    # path issues when renaming parent directories
    all_paths = []
    
    try:
        for root, dirs, files in os.walk(directory):
            # Collect all file paths
            for file in files:
                file_path = os.path.join(root, file)
                all_paths.append(('file', file_path, file))
                results['total_files'] += 1
            
            # Collect all directory paths (if we're renaming directories)
            if rename_dirs:
                for dir_name in dirs:
                    dir_path = os.path.join(root, dir_name)
                    all_paths.append(('dir', dir_path, dir_name))
                    results['total_dirs'] += 1
    
    except Exception as e:
        results['errors'].append(f"Error scanning directory: {e}")
        return results
    
    # Sort by depth (deepest first) to handle directory renames properly
    all_paths.sort(key=lambda x: x[1].count(os.sep), reverse=True)
    
    # Process each path
    for path_type, full_path, name in all_paths:
        if not os.path.exists(full_path):
            # Path might have been moved due to parent directory rename
            continue
            
        has_non_ascii = has_non_ascii_chars(name)
        
        if path_type == 'file':
            if has_non_ascii:
                results['files_with_non_ascii'] += 1
                
                # Convert filename to ASCII
                ascii_filename = convert_filename_to_ascii(name, **conversion_options)
                
                # Attempt to rename
                rename_result = rename_file(full_path, ascii_filename, dry_run)
                results['files_processed'].append(rename_result)
                
                if rename_result['success'] and rename_result['was_renamed']:
                    results['files_renamed'] += 1
                elif rename_result['error']:
                    results['errors'].append(f"Failed to rename {full_path}: {rename_result['error']}")
        
        elif path_type == 'dir' and rename_dirs:
            if has_non_ascii:
                results['dirs_with_non_ascii'] += 1
                
                # Convert directory name to ASCII
                ascii_dirname = convert_filename_to_ascii(name, **conversion_options)
                
                # Attempt to rename
                rename_result = rename_file(full_path, ascii_dirname, dry_run)
                results['dirs_processed'].append(rename_result)
                
                if rename_result['success'] and rename_result['was_renamed']:
                    results['dirs_renamed'] += 1
                elif rename_result['error']:
                    results['errors'].append(f"Failed to rename {full_path}: {rename_result['error']}")
    
    return results


def print_results(results, dry_run=False):
    """Print the results of the rename operation."""
    action_word = "Would rename" if dry_run else "Renamed"
    
    print(f"\nSummary:")
    print(f"Total files scanned: {results['total_files']}")
    if results['total_dirs'] > 0:
        print(f"Total directories scanned: {results['total_dirs']}")
    
    print(f"Files with non-ASCII names: {results['files_with_non_ascii']}")
    if results['dirs_with_non_ascii'] > 0:
        print(f"Directories with non-ASCII names: {results['dirs_with_non_ascii']}")
    
    if results['files_with_non_ascii'] == 0 and results['dirs_with_non_ascii'] == 0:
        print("✓ All files and directories already have ASCII-only names!")
        return
    
    print(f"{action_word} files: {results['files_renamed']}")
    if results['dirs_renamed'] > 0:
        print(f"{action_word} directories: {results['dirs_renamed']}")
    
    # Show renamed files
    if results['files_processed']:
        print(f"\nFiles {action_word.lower()}:")
        print("-" * 60)
        for result in results['files_processed']:
            if result['was_renamed']:
                print(f"  {result['old_filename']} → {result['new_filename']}")
                print(f"    {result['old_path']}")
    
    # Show renamed directories
    if results['dirs_processed']:
        print(f"\nDirectories {action_word.lower()}:")
        print("-" * 60)
        for result in results['dirs_processed']:
            if result['was_renamed']:
                print(f"  {result['old_filename']} → {result['new_filename']}")
                print(f"    {result['old_path']}")
    
    # Show errors
    if results['errors']:
        print(f"\nErrors:")
        print("-" * 30)
        for error in results['errors']:
            print(f"  ⚠ {error}")


def main():
    parser = argparse.ArgumentParser(
        description="Recursively rename files with non-ASCII characters to ASCII equivalents"
    )
    parser.add_argument(
        "directory", 
        nargs="?", 
        default=".", 
        help="Directory to search (default: current directory)"
    )
    parser.add_argument(
        "--include-dirs", 
        action="store_true", 
        help="Also rename directories with non-ASCII names"
    )
    parser.add_argument(
        "--dry-run", 
        action="store_true", 
        help="Show what would be renamed without making changes"
    )
    parser.add_argument(
        "--no-custom-map", 
        action="store_true", 
        help="Don't use custom Unicode-to-ASCII mapping"
    )
    parser.add_argument(
        "--no-accent-removal", 
        action="store_true", 
        help="Don't remove accents from characters"
    )
    parser.add_argument(
        "--fallback-char", 
        default="_", 
        help="Character to use for unconvertible Unicode (default: _)"
    )
    parser.add_argument(
        "--no-clean-problematic", 
        action="store_true", 
        help="Don't clean characters that are problematic in filenames"
    )
    
    args = parser.parse_args()
    
    # Validate directory
    if not os.path.isdir(args.directory):
        print(f"Error: '{args.directory}' is not a valid directory", file=sys.stderr)
        sys.exit(1)
    
    print(f"Scanning directory: {os.path.abspath(args.directory)}")
    if args.dry_run:
        print("DRY RUN - No files will be actually renamed")
    print("=" * 60)
    
    # Prepare conversion options
    conversion_options = {
        'use_custom_map': not args.no_custom_map,
        'remove_accents_flag': not args.no_accent_removal,
        'fallback_char': args.fallback_char,
        'clean_problematic': not args.no_clean_problematic
    }
    
    # Find and rename files
    results = find_and_rename_files(
        args.directory,
        rename_dirs=args.include_dirs,
        dry_run=args.dry_run,
        **conversion_options
    )
    
    # Print results
    print_results(results, args.dry_run)
    
    # Exit with appropriate code
    if results['errors']:
        sys.exit(2)  # Errors occurred
    elif results['files_with_non_ascii'] > 0 or results['dirs_with_non_ascii'] > 0:
        sys.exit(1 if not args.dry_run else 0)  # Non-ASCII files found
    else:
        sys.exit(0)  # All clean


if __name__ == "__main__":
    main()
