# Contributing

## Developer dependencies

Install developer tools:

    pip install -r requirements_dev.txt

## Testing

Run tests in isolated environments:

    tox -e py39
    tox

Tun pytest directly in your active environment:

    pytest

## Code Style

Code conforms to the `black` and PEP8 style guides. Before checking in code, please run the linters:

    black .
    flake8

These are tested by the 'lint' tox environment:

    tox -e lint
