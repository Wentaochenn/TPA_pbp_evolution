from pathlib import Path

def command(executable, posterior, target, output):
    return [executable, '-lowMem', 'true', '-burnin', '10', '-height', 'CA',
            '-topology', 'target', '-target', target, posterior, output]

def require_tree(path, log_path):
    path, log_path = Path(path), Path(log_path)
    if not path.is_file() or not path.stat().st_size:
        tail = '\n'.join(log_path.read_text(errors='replace').splitlines()[-55:]) if log_path.exists() else '<missing log>'
        raise ValueError(f'TreeAnnotator returned without a tree. Log: {log_path}\n{tail}')
