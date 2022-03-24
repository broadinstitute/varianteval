import varianteval.engine.config_utils as config_utils
import argparse
from varianteval.core.io import vcf

# ------ CLI ------
parser = argparse.ArgumentParser(description='Generate an SV callset report')
parser.add_argument('--config', help='Report YAML config')
args = parser.parse_args()
# -----------------

config = config_utils.load_config(args.config, config_type=config_utils.CONFIG_TYPE.REPORT)
callset = vcf.vcf2callset(config.vcf)
callset.report()
