# 公開配布版のコード構成

起動方法と現在値はルートREADME.mdを参照してください。

- portable_launcher.py: portable entry and worker dispatch.
- viewer/xray_viewer.py: input, projection, settings and recent simulation speed.
- viewer/native_engine.py: calculation worker.
- viewer/native_core.py: ctypes interface to runtime/guidewire_core.dll.
- native/bridge.cpp: C interface.
- native/kernel.h: insertion, beam force, system assembly and free motion.
- native/scene.h: collision, friction constraints and correction.
- native/generated/: extracted SOFA/BeamAdapter function bodies and manifests.
- upstream/: corresponding original source, math headers and licenses.

ルートbuild_native.ps1で計算核、build_portable.ps1でexeを再ビルドします。
ソース版だけの場合は計算核のビルドとPython依存のインストールが必要です。
起動用ZIPにはビルド済み計算核とPython実行環境を含みます。