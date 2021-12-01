from pathlib import Path

def file_info(path):
    if not path.exists():
        return None
    return {
        'name': path.name,
        'size': path.stat().st_size
    }

def collect_files(dir, basename):
    return {
        'TSV': file_info(dir / f'{basename}.tsv'),
        'PARQUET': file_info(dir / f'{basename}.parquet')
    }

def parse_dl_files(directory):
    result = {}
    dir = Path(directory)
    result['annot'] = collect_files(dir, 'deg_analyses')
    result['all'] = collect_files(dir, 'all_viruses')
    result['virus'] = {}
    for tsv in sorted(dir.glob('*.tsv')):
        if tsv.stem in set(['all_viruses', 'deg_analyses']):
            continue
        virus = tsv.stem.replace('_', ' ')
        result['virus'][virus] = collect_files(dir, tsv.stem)
    return result