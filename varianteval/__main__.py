import logging

import click
import click_log

import sys
from datetime import datetime

import varianteval
from . import __version__


logger = logging.getLogger(__name__)
click_log.basic_config(logger)
logger.handlers[0].formatter = logging.Formatter(
    "[%(levelname)s %(asctime)s %(name)8s] %(message)s", "%Y-%m-%d %H:%M:%S"
)


@click.group(name="varianteval")
def main_entry():
    logger.info("Invoked via: varianteval %s", " ".join(sys.argv))


@main_entry.command()
@click_log.simple_verbosity_option(logger)
def version():
    """Print the version of varianteval."""
    click.echo(__version__)


if __name__ == "__main__":
    main_entry()