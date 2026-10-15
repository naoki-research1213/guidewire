"""このリポジトリの原ソースから再抽出する。vendorやSOFAインストール不要。"""
from pathlib import Path
import subprocess,sys
HERE=Path(__file__).resolve().parent
for name in ('native','solver','contact_solver','collision','bvh','mesh','mapping','barycentric','contact_model','history','identifier','subset'):
    subprocess.run([sys.executable,str(HERE/f'extract_{name}.py')],check=True)
print('Regenerated from tracked upstream sources')
