#!/usr/bin/env python3
"""
Script to replace Unicode characters in text files with their ASCII equivalents.

There are various options, but *.txt will fix all of the txt files in a directory
This routine is not recursive in terms of directory structure
"""

import os
import sys
import argparse
import unicodedata
import re
from pathlib import Path


# Custom mapping for common Unicode characters to ASCII equivalents
UNICODE_TO_ASCII_MAP = {
    # Quotation marks
    '\u2018': "'",      # LEFT SINGLE QUOTATION MARK
    '\u2019': "'",      # RIGHT SINGLE QUOTATION MARK
    '\u201A': "'",      # SINGLE LOW-9 QUOTATION MARK
    '\u201B': "'",      # SINGLE HIGH-REVERSED-9 QUOTATION MARK
    '\u201C': '"',      # LEFT DOUBLE QUOTATION MARK
    '\u201D': '"',      # RIGHT DOUBLE QUOTATION MARK
    '\u201E': '"',      # DOUBLE LOW-9 QUOTATION MARK
    '\u201F': '"',      # DOUBLE HIGH-REVERSED-9 QUOTATION MARK
    '\u2039': '<',      # SINGLE LEFT-POINTING ANGLE QUOTATION MARK
    '\u203A': '>',      # SINGLE RIGHT-POINTING ANGLE QUOTATION MARK
    '\u00AB': '<<',     # LEFT-POINTING DOUBLE ANGLE QUOTATION MARK
    '\u00BB': '>>',     # RIGHT-POINTING DOUBLE ANGLE QUOTATION MARK
    
    # Dashes and hyphens
    '\u2010': '-',      # HYPHEN
    '\u2011': '-',      # NON-BREAKING HYPHEN
    '\u2012': '-',      # FIGURE DASH
    '\u2013': '-',      # EN DASH
    '\u2014': '--',     # EM DASH
    '\u2015': '--',     # HORIZONTAL BAR
    '\u2212': '-',      # MINUS SIGN
    
    # Spaces and breaks
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
    
    # Dots and bullets
    '\u2022': '*',      # BULLET
    '\u2023': '>',      # TRIANGULAR BULLET
    '\u2024': '.',      # ONE DOT LEADER
    '\u2025': '..',     # TWO DOT LEADER
    '\u2026': '...',    # HORIZONTAL ELLIPSIS
    '\u00B7': '*',      # MIDDLE DOT
    '\u2219': '*',      # BULLET OPERATOR
    
    # Mathematical symbols
    '\u00B1': '+-',     # PLUS-MINUS SIGN
    '\u00D7': 'x',      # MULTIPLICATION SIGN
    '\u00F7': '/',      # DIVISION SIGN
    '\u2260': '!=',     # NOT EQUAL TO
    '\u2264': '<=',     # LESS-THAN OR EQUAL TO
    '\u2265': '>=',     # GREATER-THAN OR EQUAL TO
    '\u00B0': 'deg',    # DEGREE SIGN
    '\u2032': "'",      # PRIME
    '\u2033': '"',      # DOUBLE PRIME
    '\u221E': 'inf',    # INFINITY
    
    # Currency symbols
    '\u00A2': 'c',      # CENT SIGN
    '\u00A3': 'GBP',    # POUND SIGN
    '\u00A4': '$',      # GENERIC CURRENCY SYMBOL
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


def normalize_unicode(text, method='nfd'):
    """
    Normalize Unicode text using different methods.
    
    Args:
        text (str): Input text
        method (str): Normalization method ('nfd', 'nfc', 'nfkd', 'nfkc')
    
    Returns:
        str: Normalized text
    """
    if method.upper() == 'NFD':
        return unicodedata.normalize('NFD', text)
    elif method.upper() == 'NFC':
        return unicodedata.normalize('NFC', text)
    elif method.upper() == 'NFKD':
        return unicodedata.normalize('NFKD', text)
    elif method.upper() == 'NFKC':
        return unicodedata.normalize('NFKC', text)
    else:
        return text


def remove_accents(text):
    """Remove accents from characters using Unicode normalization."""
    # Normalize to NFD (decomposed form) and filter out combining characters
    normalized = unicodedata.normalize('NFD', text)
    ascii_text = ''.join(
        char for char in normalized 
        if unicodedata.category(char) != 'Mn'  # Mn = Nonspacing_Mark (accents)
    )
    return ascii_text


def replace_unicode_with_ascii(text, 
                              use_custom_map=True, 
                              remove_accents_flag=True,
                              transliterate=True,
                              fallback_char='?',
                              preserve_newlines=True):
    """
    Replace Unicode characters with ASCII equivalents using multiple strategies.
    
    Args:
        text (str): Input text
        use_custom_map (bool): Use custom Unicode-to-ASCII mapping
        remove_accents_flag (bool): Remove accents from characters
        transliterate (bool): Try to transliterate remaining characters
        fallback_char (str): Character to use when no conversion is possible
        preserve_newlines (bool): Keep newline characters as-is
    
    Returns:
        tuple: (converted_text, conversion_stats)
    """
    original_length = len(text)
    conversions = 0
    unconvertible = 0
    
    result = text
    
    # Step 1: Apply custom mapping
    if use_custom_map:
        for unicode_char, ascii_equiv in UNICODE_TO_ASCII_MAP.items():
            if unicode_char in result:
                count = result.count(unicode_char)
                result = result.replace(unicode_char, ascii_equiv)
                conversions += count
    
    # Step 2: Remove accents (decompose and remove combining marks)
    if remove_accents_flag:
        before_accent_removal = result
        result = remove_accents(result)
        # Count characters that were actually changed
        for i, (old_char, new_char) in enumerate(zip(before_accent_removal, result)):
            if old_char != new_char and ord(old_char) > 127:
                conversions += 1
    
    # Step 3: Handle remaining non-ASCII characters
    final_result = []
    for char in result:
        if ord(char) <= 127:
            # ASCII character, keep as-is
            final_result.append(char)
        elif preserve_newlines and char in '\n\r\t':
            # Preserve common whitespace
            final_result.append(char)
        elif transliterate:
            # Try to find ASCII representation
            try:
                # Try NFKD normalization first
                normalized = unicodedata.normalize('NFKD', char)
                ascii_version = ''.join(
                    c for c in normalized 
                    if ord(c) <= 127 or unicodedata.category(c) != 'Mn'
                )
                if ascii_version and all(ord(c) <= 127 for c in ascii_version):
                    final_result.append(ascii_version)
                    conversions += 1
                else:
                    # Use fallback character
                    final_result.append(fallback_char)
                    unconvertible += 1
            except:
                final_result.append(fallback_char)
                unconvertible += 1
        else:
            # Use fallback character
            final_result.append(fallback_char)
            unconvertible += 1
    
    final_text = ''.join(final_result)
    
    stats = {
        'original_length': original_length,
        'final_length': len(final_text),
        'conversions': conversions,
        'unconvertible': unconvertible,
        'ascii_percentage': sum(1 for c in final_text if ord(c) <= 127) / len(final_text) * 100 if final_text else 100
    }
    
    return final_text, stats


def has_non_ascii_chars(text):
    """Check if text contains any non-ASCII characters."""
    return any(ord(char) > 127 for char in text)


def create_backup(input_path, backup_dir='original'):
    """
    Create a backup copy of the input file.
    
    Args:
        input_path (str): Path to the file to backup
        backup_dir (str): Directory to store backups
    
    Returns:
        str: Path to backup file, or None if backup failed
    """
    try:
        os.makedirs(backup_dir, exist_ok=True)
        filename = os.path.basename(input_path)
        backup_path = os.path.join(backup_dir, filename)
        
        # If backup already exists, add a number suffix
        counter = 1
        base_backup_path = backup_path
        while os.path.exists(backup_path):
            name, ext = os.path.splitext(base_backup_path)
            backup_path = f"{name}_{counter}{ext}"
            counter += 1
        
        # Copy the file
        import shutil
        shutil.copy2(input_path, backup_path)
        return backup_path
        
    except Exception as e:
        return None


def convert_file(input_path, output_path=None, encoding='utf-8', create_backup_flag=True, backup_dir='original', **conversion_options):
    """
    Convert a single file from Unicode to ASCII.
    
    Args:
        input_path (str): Path to input file
        output_path (str): Path to output file (None for in-place)
        encoding (str): File encoding
        create_backup_flag (bool): Whether to create backup copies
        backup_dir (str): Directory for backup files
        **conversion_options: Options for Unicode conversion
    
    Returns:
        dict: Conversion results
    """
    results = {
        'input_path': input_path,
        'output_path': output_path or input_path,
        'success': False,
        'error': None,
        'stats': None,
        'backup_path': None,
        'had_non_ascii': False,
        'backup_created': False
    }
    
    try:
        # Read input file
        with open(input_path, 'r', encoding=encoding, errors='replace') as f:
            original_text = f.read()
        
        # Check if file has non-ASCII characters
        results['had_non_ascii'] = has_non_ascii_chars(original_text)
        
        # Convert text
        converted_text, stats = replace_unicode_with_ascii(original_text, **conversion_options)
        
        # Create backup if file has non-ASCII characters and backup is requested
        if results['had_non_ascii'] and create_backup_flag:
            backup_path = create_backup(input_path, backup_dir)
            if backup_path:
                results['backup_path'] = backup_path
                results['backup_created'] = True
        
        # Only write output if there were actual changes or if forced
        if results['had_non_ascii'] or stats['conversions'] > 0:
            output_file = output_path or input_path
            with open(output_file, 'w', encoding='ascii', errors='replace') as f:
                f.write(converted_text)
        
        results['success'] = True
        results['stats'] = stats
        
    except Exception as e:
        results['error'] = str(e)
    
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Convert Unicode characters in text files to ASCII equivalents"
    )
    parser.add_argument(
        "files", 
        nargs="+", 
        help="Input file(s) to convert"
    )
    parser.add_argument(
        "--output-dir", 
        help="Output directory (default: overwrite input files)"
    )
    parser.add_argument(
        "--suffix", 
        default="", 
        help="Suffix to add to output filenames (e.g., '_ascii')"
    )
    parser.add_argument(
        "--encoding", 
        default="utf-8", 
        help="Input file encoding (default: utf-8)"
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
        "--no-transliterate", 
        action="store_true", 
        help="Don't attempt transliteration"
    )
    parser.add_argument(
        "--fallback-char", 
        default="?", 
        help="Character to use for unconvertible Unicode (default: ?)"
    )
    parser.add_argument(
        "--dry-run", 
        action="store_true", 
        help="Show what would be converted without making changes"
    )
    parser.add_argument(
        "--no-backup", 
        action="store_true", 
        help="Don't create backup copies of files with non-ASCII characters"
    )
    parser.add_argument(
        "--backup-dir", 
        default="original", 
        help="Directory for backup files (default: original)"
    )
    
    args = parser.parse_args()
    
    # Prepare conversion options
    conversion_options = {
        'use_custom_map': not args.no_custom_map,
        'remove_accents_flag': not args.no_accent_removal,
        'transliterate': not args.no_transliterate,
        'fallback_char': args.fallback_char,
        'preserve_newlines': True
    }
    
    total_files = len(args.files)
    successful_conversions = 0
    
    for file_path in args.files:
        print(f"Processing: {file_path}")
        
        if not os.path.isfile(file_path):
            print(f"  Error: File not found", file=sys.stderr)
            continue
        
        # Determine output path
        if args.output_dir:
            os.makedirs(args.output_dir, exist_ok=True)
            filename = os.path.basename(file_path)
            if args.suffix:
                name, ext = os.path.splitext(filename)
                filename = f"{name}{args.suffix}{ext}"
            output_path = os.path.join(args.output_dir, filename)
        elif args.suffix:
            name, ext = os.path.splitext(file_path)
            output_path = f"{name}{args.suffix}{ext}"
        else:
            output_path = None  # In-place conversion
        
        if args.dry_run:
            # Just analyze what would be converted
            try:
                with open(file_path, 'r', encoding=args.encoding, errors='replace') as f:
                    original_text = f.read()
                
                has_non_ascii = has_non_ascii_chars(original_text)
                converted_text, stats = replace_unicode_with_ascii(original_text, **conversion_options)
                
                if not has_non_ascii:
                    print(f"  ✓ File contains only ASCII characters - no conversion needed")
                else:
                    print(f"  Would convert: {stats['conversions']} characters")
                    print(f"  Unconvertible: {stats['unconvertible']} characters")
                    print(f"  Final ASCII %: {stats['ascii_percentage']:.1f}%")
                    if not args.no_backup:
                        print(f"  Would create backup in: {args.backup_dir}/")
                    if output_path:
                        print(f"  Would save to: {output_path}")
                    else:
                        print(f"  Would overwrite: {file_path}")
            except Exception as e:
                print(f"  Error during analysis: {e}", file=sys.stderr)
        else:
            # Actually convert the file
            results = convert_file(
                file_path, 
                output_path, 
                args.encoding, 
                create_backup_flag=not args.no_backup,
                backup_dir=args.backup_dir,
                **conversion_options
            )
            
            if results['success']:
                if not results['had_non_ascii']:
                    print(f"  ✓ File contains only ASCII characters - no conversion needed")
                else:
                    stats = results['stats']
                    print(f"  ✓ Converted: {stats['conversions']} characters")
                    print(f"  ✓ Unconvertible: {stats['unconvertible']} characters")
                    print(f"  ✓ Final ASCII %: {stats['ascii_percentage']:.1f}%")
                    if results['backup_created']:
                        print(f"  ✓ Backup created: {results['backup_path']}")
                    elif not args.no_backup:
                        print(f"  ! Backup creation failed")
                    print(f"  ✓ Saved to: {results['output_path']}")
                successful_conversions += 1
            else:
                print(f"  ✗ Error: {results['error']}", file=sys.stderr)
        
        print()
    
    if not args.dry_run:
        print(f"Successfully converted {successful_conversions}/{total_files} files")
    
    # Exit with error code if not all files were processed successfully
    sys.exit(0 if successful_conversions == total_files else 1)


if __name__ == "__main__":
    main()
