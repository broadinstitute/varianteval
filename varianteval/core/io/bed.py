import functools
import bisect
from varianteval.core.constants import *
from varianteval.core.variant.intervals import GenomeInterval

class BedRecord:
    def __init__(self, interval_a, interval_b=None, sv_type=None, qual=None, gt=None):
        self.format = BED_TYPE.BEDPE if interval_b is not None else BED_TYPE.BED
        self.interval_a = interval_a
        self.interval_b = interval_b
        self.sv_type = sv_type
        self.qual = qual
        self.gt = gt

    @staticmethod
    def parse_bed_line(line, bed_file_type=BED_TYPE.BEDPE):
        """
        Parses a single line of a BED (one interval) or BEDPE (two intervals) file
        Optional auxiliary (i.e. non-positional) fields are assumed to be given in the following order:
            SV type, SV quality score, genotype
        """
        fields = line.strip().split()
        assert len(fields) >= 3, "Unexpected number of fields in BED: %s" % line
        chr_name, start, end = fields[:3]
        bed_record = BedRecord(GenomeInterval(chr_name, int(start), int(end)))
        if bed_file_type == BED_TYPE.BEDPE:
            assert len(fields) >= 6, "Unexpected number of fields in BEDPE: %s" % line
            chr_b, start_b, end_b = fields[3:6]
            bed_record.interval_b = GenomeInterval(chr_b, int(start_b), int(end_b))
        req_fields = 3 if bed_file_type == BED_TYPE.BED else 6
        bed_record.sv_type = SV_TYPE_ENCODING[fields[req_fields] if len(fields) > req_fields else 'UNK']
        bed_record.qual = fields[req_fields + 1] if len(fields) > req_fields + 1 else None
        bed_record.gt = BED2GT_ENCODING[fields[req_fields + 2]] if len(fields) > req_fields + 2 else GT.UNK
        # TODO: strand vs zygosity
        return bed_record

    def get_sv_type(self, to_vcf_format=False):
        if to_vcf_format:
            return "<%s>" % self.sv_type.name
        return self.sv_type.name

    def get_qual(self, convert_none=False):
        if self.qual is None:
            return 0 if convert_none else None
        return self.qual

    def get_gt(self):
        return self.gt

    def get_sv_type_with_gt(self):
        return "%s-%s" % (self.sv_type.name, self.get_gt().name)

    def __str__(self):
        return "%s: %s, %s" % (self.sv_type, str(self.interval_a), str(self.interval_b))

    def get_name(self):
        return "%s_%s_%s" % (self.sv_type, str(self.interval_a), str(self.interval_b))

    def to_bedpe(self):
        assert self.format == BED_TYPE.BEDPE
        return "%s\t%s\t%s\t%s\t%s\t%s" % (self.interval_a.chr_name, self.interval_a.start, self.interval_a.end,
                                           self.interval_b.chr_name, self.interval_b.start, self.interval_b.end)

    def to_bedpe_aux(self):
        assert self.format == BED_TYPE.BEDPE
        return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (self.interval_a.chr_name,
                                                       self.interval_a.start, self.interval_a.end,
                                                       self.interval_b.chr_name,
                                                       self.interval_b.start, self.interval_b.end,
                                                       self.sv_type, self.qual,
                                                       GT2BED_ENCODING[self.gt])

    def to_bed(self):
        assert self.format == BED_TYPE.BED
        return "%s\t%s\t%s" % (self.interval_a.chr_name, self.interval_a.start, self.interval_a.end)

    @staticmethod
    def get_bedpe_header():
        return '#chrom1\tstart1\tstop1\tchrom2\tstart2\tstop2'

    @staticmethod
    def get_bedpe_aux_header():
        return '#chrom1\tstart1\tstop1\tchrom2\tstart2\tstop2\ttype\tqual\tgt'

    @staticmethod
    def compare(rec1, rec2):
        return rec1.interval_a.__lt__(rec2.interval_a)

    @staticmethod
    def compare_by_qual(rec1, rec2):
        return rec1.get_qual(convert_none=True) - rec2.get_qual(convert_none=True)

    def __lt__(self, rec):
        return self.interval_a.__lt__(rec.interval_a)

def bed_iter(bed_fname, bed_file_type=BED_TYPE.BEDPE, keep_chrs=None, keep_sv_types=None):
    def skip_record(bed_record, chr_whitelist=None, type_whitelist=None):
        return (type_whitelist is not None and bed_record.sv_type not in type_whitelist) or \
               (chr_whitelist is not None and (bed_record.interval_a.chr_name not in chr_whitelist or
                                               bed_record.interval_b is not None and
                                               bed_record.interval_b.chr_name not in chr_whitelist))

    with open(bed_fname, 'r') as bed_file:
        for line in bed_file:
            if line.startswith('#') or line.isspace():
                continue
            record = BedRecord.parse_bed_line(line, bed_file_type)
            if not skip_record(record, keep_chrs, keep_sv_types):
                yield record

def load_bed(bed_fname, bed_file_type=BED_TYPE.BEDPE, sort=True, keep_chrs=None, keep_sv_types=None):
    records = []
    for record in bed_iter(bed_fname, bed_file_type, keep_chrs=keep_chrs, keep_sv_types=keep_sv_types):
        records.append(record)
    if sort:
        records = sorted(records, key=functools.cmp_to_key(BedRecord.compare))
    return records

def write_bedpe(out_bed_file, bed_records):
    bed_file = open(out_bed_file, 'w')
    bed_file.write(BedRecord.get_bedpe_aux_header() + "\n")
    for rec in bed_records:
        bed_file.write(rec.to_bedpe_aux() + "\n")
    bed_file.close()
