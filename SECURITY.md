# Security Policy

## Reporting a Vulnerability

If you discover a security vulnerability in BBCN, please **do not** open a public issue.

Instead, please email the maintainer with:
- Description of the vulnerability
- Steps to reproduce (if applicable)
- Potential impact
- Suggested fix (if you have one)

We take security seriously and will work with you to understand and address the issue promptly.

## Known Limitations

- BBCN is a **research tool**, not a clinical decision-making system
- Model predictions should not be used alone for patient treatment decisions
- All outputs should be validated independently
- See `docs/SCOPE.md` for detailed limitations and disclaimers

## Dependencies

This project depends on:
- `numpy` and `pandas` for computation
- `matplotlib` for figure generation (optional)

We monitor these dependencies for security issues. Please report if you identify vulnerabilities in any dependency.