# Contributing to BBCN

Thank you for your interest in contributing to the Boolean Breast-Cancer Network project! This document outlines how to report issues, suggest improvements, and contribute code.

## Reporting Issues

### Bug Reports
If you find a bug, please open an issue with:
- **Title**: Brief, descriptive summary
- **Environment**: Python version, OS, any relevant versions
- **Reproduction steps**: Exact commands to reproduce the issue
- **Expected vs. actual behavior**: What should happen vs. what does
- **Error output**: Full traceback if applicable

### Feature Requests
Suggest new features by opening an issue with:
- **Title**: Clear description of the feature
- **Motivation**: Why this feature would be useful
- **Proposed implementation** (if you have ideas)

## Contributing Code

### Setup
1. Fork the repository
2. Clone your fork: `git clone https://github.com/YOUR_USERNAME/Systems_Biology_Journal_2026.git`
3. Create a feature branch: `git checkout -b feature/your-feature-name`
4. Install dependencies: `pip install -r requirements.txt`

### Code Guidelines
- **Python style**: Follow PEP 8 conventions
- **Comments**: Include docstrings for all functions and classes
- **Testing**: Ensure `python reproduce_all.py` runs without errors
- **Commit messages**: Write clear, concise messages (e.g., "Fix typo in controller.py" or "Add validation to seed_from_signatures.py")

### Submitting Changes
1. Push your branch to your fork
2. Open a pull request with:
   - **Title**: Clear description of changes
   - **Description**: Explain what changed and why
   - **Verification**: Confirm that `python reproduce_all.py` completes successfully

## Code of Conduct
Be respectful and constructive in all interactions. We're building scientific tools together.

## Questions?
Open an issue with the `question` label or check the existing docs in `docs/`.

Thank you for contributing!