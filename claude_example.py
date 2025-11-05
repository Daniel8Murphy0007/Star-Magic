#!/usr/bin/env python3
"""
Claude Integration Example for Star-Magic Project

This script demonstrates how to use Anthropic's Claude API to analyze
the Star-Magic theoretical physics content.

Requirements:
    - anthropic package (install via: pip install -r requirements.txt)
    - ANTHROPIC_API_KEY environment variable set

Usage:
    export ANTHROPIC_API_KEY='your-api-key-here'
    python claude_example.py
"""

import os
from anthropic import Anthropic


def main():
    """
    Example function demonstrating Claude API usage with Star-Magic content.
    """
    # Initialize the Anthropic client
    # The API key should be set in the ANTHROPIC_API_KEY environment variable
    api_key = os.environ.get("ANTHROPIC_API_KEY")
    
    if not api_key:
        print("Error: ANTHROPIC_API_KEY environment variable not set.")
        print("Please set it with: export ANTHROPIC_API_KEY='your-api-key-here'")
        return
    
    client = Anthropic(api_key=api_key)
    
    # Example: Ask Claude about the Star-Magic theory
    print("Star-Magic Claude Integration Example")
    print("=" * 50)
    print("\nAsking Claude about the Unified Quantum Field Equation...\n")
    
    try:
        # Create a message asking about the Star-Magic theory
        message = client.messages.create(
            model="claude-3-5-sonnet-20241022",
            max_tokens=1024,
            messages=[
                {
                    "role": "user",
                    "content": "Explain what a unified field equation aims to achieve in physics."
                }
            ]
        )
        
        # Print the response
        print("Claude's Response:")
        print("-" * 50)
        print(message.content[0].text)
        print("-" * 50)
        print("\n✓ Claude integration successful!")
        print("\nYou can now use this client to analyze Star-Magic content,")
        print("ask questions about the theory, or process the theoretical framework.")
        
    except Exception as e:
        print(f"Error communicating with Claude API: {e}")
        print("\nMake sure your API key is valid and you have internet connectivity.")


if __name__ == "__main__":
    main()
