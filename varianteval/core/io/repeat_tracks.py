import gzip

def repeat_track_iter(repeat_track_fname):
    """ Iterator for a RepeatMasker track file (assumes .txt.gz format) """
    with gzip.open(repeat_track_fname, "r") as track:
        for line in track:
            entries = line.strip().split()
            chr = entries[0]
            start = int(entries[1])
            end = int(entries[2])
            repeat_name = entries[3]
            yield chr, start, end, repeat_name
