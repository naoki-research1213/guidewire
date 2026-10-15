"""抽出入力はGitに保存された原ソース。外部SOFAのインストールは不要。"""
from pathlib import Path
import shutil
ROOT=Path(__file__).resolve().parents[1]
def source_file(path):
    path=Path(path)
    if 'vendor' not in path.parts:return path
    found=list((ROOT/'upstream').rglob(path.name))
    if len(found)!=1:raise FileNotFoundError(f'Expected one tracked source for {path.name}, found {found}')
    return found[0]
def copy_source(path,destination):
    path=source_file(path);destination=Path(destination)
    if path.resolve()!=destination.resolve():shutil.copy2(path,destination)
