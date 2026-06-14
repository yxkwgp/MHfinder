# Contributing

Thank you for improving MHfinder.

## Development setup

```bash
git clone https://github.com/Fun-Gene/MHfinder.git
cd MHfinder
pip install -e .[analysis,dev]
python -m unittest discover -s tests
```

PLINK and R are required only for end-to-end genomic workflows, not for lightweight unit tests.

## Reporting issues

Please include:

- MHfinder version or Git commit,
- operating system,
- Python/R/PLINK versions,
- exact command run,
- relevant input-file format summary, and
- full error message.

Do not upload restricted genotype data unless it is public and explicitly shareable.
