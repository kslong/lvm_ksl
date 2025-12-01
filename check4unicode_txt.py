#!/usr/bin/env python3
"""
Script to check text files for non-ASCII characters and identify what they are.

This is a belt and suspenders routine; unicode2ascii.py fixes things
"""

import os
import sys
import argparse
import unicodedata
from collections import defaultdict


def get_char_info(char):
    """Get detailed information about a Unicode character."""
    code_point = ord(char)
    try:
        name = unicodedata.name(char)
    except ValueError:
        name = "UNNAMED"
    
    category = unicodedata.category(char)
    return {
        'char': char,
        'code_point': code_point,
        'hex': f'U+{code_point:04X}',
        'name': name,
        'category': category,
        'bytes_utf8': char.encode('utf-8')
    }


def analyze_text_file(file_path, encoding='utf-8'):
    """
    Analyze a text file for non-ASCII characters.
    
    Args:
        file_path (str): Path to the text file
        encoding (str): File encoding to use
    
    Returns:
        dict: Analysis results
    """
    results = {
        'file_path': file_path,
        'total_chars': 0,
        'non_ascii_chars': 0,
        'non_ascii_positions': [],
        'char_frequency': defaultdict(int),
        'lines_with_non_ascii': set(),
        'encoding_used': encoding,
        'error': None
    }
    
    try:
        with open(file_path, 'r', encoding=encoding, errors='replace') as f:
            line_num = 1
            char_pos = 0
            
            for line in f:
                col_num = 1
                for char in line:
                    results['total_chars'] += 1
                    
                    if ord(char) > 127:  # Non-ASCII character
                        results['non_ascii_chars'] += 1
                        results['char_frequency'][char] += 1
                        results['lines_with_non_ascii'].add(line_num)
                        
                        char_info = get_char_info(char)
                        char_info.update({
                            'line': line_num,
                            'column': col_num,
                            'position': char_pos
                        })
                        results['non_ascii_positions'].append(char_info)
                    
                    char_pos += 1
                    col_num += 1
                
                line_num += 1
                
    except UnicodeDecodeError as e:
        results['error'] = f"Encoding error: {e}"
    except FileNotFoundError:
        results['error'] = f"File not found: {file_path}"
    except PermissionError:
        results['error'] = f"Permission denied: {file_path}"
    except Exception as e:
        results['error'] = f"Unexpected error: {e}"
    
    return results


def format_char_display(char_info):
    """Format character information for display."""
    char = char_info['char']
    
    # Handle special cases for display
    if char == '\n':
        display_char = '\\n'
    elif char == '\t':
        display_char = '\\t'
    elif char == '\r':
        display_char = '\\r'
    elif unicodedata.category(char) in ['Cc', 'Cf', 'Cs', 'Co', 'Cn']:
        # Control characters, format characters, etc.
        display_char = f'[{char_info["hex"]}]'
    else:
        display_char = char
    
    return display_char


def print_analysis(results, show_positions=True, max_positions=50):
    """Print the analysis results in a formatted way."""
    if results['error']:
        print(f"Error analyzing file: {results['error']}", file=sys.stderr)
        return
    
    print(f"File: {results['file_path']}")
    print(f"Encoding: {results['encoding_used']}")
    print(f"Total characters: {results['total_chars']:,}")
    print(f"Non-ASCII characters: {results['non_ascii_chars']:,}")
    
    if results['non_ascii_chars'] == 0:
        print("✓ File contains only ASCII characters")
        return
    
    print(f"Lines with non-ASCII characters: {len(results['lines_with_non_ascii'])}")
    print(f"Percentage non-ASCII: {results['non_ascii_chars']/results['total_chars']*100:.2f}%")
    print()
    
    # Show character frequency
    print("Non-ASCII character frequency:")
    print("-" * 40)
    for char, count in sorted(results['char_frequency'].items(), key=lambda x: x[1], reverse=True):
        char_info = get_char_info(char)
        display_char = format_char_display(char_info)
        print(f"'{display_char}' ({char_info['hex']}) - {char_info['name']} - {count} times")
    print()
    
    # Show positions if requested
    if show_positions and results['non_ascii_positions']:
        print(f"Character positions (showing first {max_positions}):")
        print("-" * 60)
        
        for i, char_info in enumerate(results['non_ascii_positions'][:max_positions]):
            display_char = format_char_display(char_info)
            print(f"Line {char_info['line']:4d}, Col {char_info['column']:3d}: "
                  f"'{display_char}' ({char_info['hex']}) - {char_info['name']}")
        
        if len(results['non_ascii_positions']) > max_positions:
            remaining = len(results['non_ascii_positions']) - max_positions
            print(f"... and {remaining} more")


def main():
    parser = argparse.ArgumentParser(
        description="Check text files for non-ASCII characters and identify them"
    )
    parser.add_argument(
        "files", 
        nargs="+", 
        help="Text file(s) to analyze"
    )
    parser.add_argument(
        "--encoding", 
        default="utf-8", 
        help="File encoding (default: utf-8)"
    )
    parser.add_argument(
        "--no-positions", 
        action="store_true", 
        help="Don't show character positions"
    )
    parser.add_argument(
        "--max-positions", 
        type=int, 
        default=50, 
        help="Maximum number of positions to show (default: 50)"
    )
    parser.add_argument(
        "--summary-only", 
        action="store_true", 
        help="Show only summary statistics"
    )
    
    args = parser.parse_args()
    
    exit_code = 0
    
    for file_path in args.files:
        if len(args.files) > 1:
            print("=" * 60)
        
        results = analyze_text_file(file_path, args.encoding)
        
        if results['error']:
            exit_code = 2
            print_analysis(results)
            continue
        
        if results['non_ascii_chars'] > 0:
            exit_code = 1
        
        if args.summary_only:
            print(f"{file_path}: {results['non_ascii_chars']} non-ASCII chars "
                  f"out of {results['total_chars']} total chars")
        else:
            print_analysis(results, 
                         show_positions=not args.no_positions,
                         max_positions=args.max_positions)
        
        if len(args.files) > 1:
            print()
    
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
