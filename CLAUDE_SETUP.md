# Claude Integration Guide

This guide explains how to use Anthropic's Claude AI with the Star-Magic theoretical physics repository.

## Installation

### Prerequisites
- Python 3.7 or higher
- pip (Python package manager)

### Install Dependencies

```bash
pip install -r requirements.txt
```

This will install the `anthropic` package and its dependencies.

## Configuration

### API Key Setup

To use Claude, you need an API key from Anthropic:

1. Sign up for an account at [console.anthropic.com](https://console.anthropic.com)
2. Generate an API key from your account settings
3. Set the API key as an environment variable:

```bash
export ANTHROPIC_API_KEY='your-api-key-here'
```

For persistent configuration, add this to your `~/.bashrc` or `~/.zshrc`:

```bash
echo 'export ANTHROPIC_API_KEY="your-api-key-here"' >> ~/.bashrc
source ~/.bashrc
```

## Usage

### Basic Example

Run the included example script:

```bash
python claude_example.py
```

This demonstrates basic Claude API usage and can be used as a template for your own integrations.

### Using Claude with Star-Magic Content

You can use Claude to:
- Analyze the Unified Quantum Field Equation
- Discuss theoretical implications
- Generate explanations of complex concepts
- Process and summarize the Star-Magic framework
- Explore mathematical derivations

Example Python code:

```python
from anthropic import Anthropic
import os

client = Anthropic(api_key=os.environ.get("ANTHROPIC_API_KEY"))

# Read content from Star Magic.md
with open("Star Magic.md", "r") as f:
    content = f.read()

# Ask Claude to analyze specific sections
message = client.messages.create(
    model="claude-3-5-sonnet-20241022",
    max_tokens=2048,
    messages=[
        {
            "role": "user",
            "content": f"Analyze this theoretical physics content: {content[:4000]}"
        }
    ]
)

print(message.content[0].text)
```

## Available Models

The anthropic package supports various Claude models:
- `claude-3-5-sonnet-20241022` (recommended for most use cases)
- `claude-3-opus-20240229` (most capable for complex reasoning)
- `claude-3-sonnet-20240229` (balanced performance)
- `claude-3-haiku-20240307` (fastest, most cost-effective)

## Resources

- [Anthropic Documentation](https://docs.anthropic.com)
- [Claude API Reference](https://docs.anthropic.com/claude/reference)
- [Star-Magic Theory Documentation](Star%20Magic.md)

## Troubleshooting

### Common Issues

**"ANTHROPIC_API_KEY environment variable not set"**
- Ensure you've exported the API key: `export ANTHROPIC_API_KEY='your-key'`
- Check if the key is set: `echo $ANTHROPIC_API_KEY`

**"API key is invalid"**
- Verify your API key at [console.anthropic.com](https://console.anthropic.com)
- Make sure there are no extra spaces or quotes in your key

**Connection errors**
- Check your internet connection
- Verify you can access anthropic.com
- Check if your firewall allows HTTPS connections

## Support

For issues related to:
- Claude API: [Anthropic Support](https://support.anthropic.com)
- Star-Magic theory: Open an issue in this repository
