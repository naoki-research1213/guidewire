# Guidewire Simulator — Windows portable snapshot

## 起動（Python・SOFAのインストール不要）

1. [起動用ZIPをダウンロード](https://github.com/naoki-research1213/guidewire/releases/latest/download/GuidewireSimulator-Windows-x64-2026-10-15.zip)。
2. ZIPを右クリックして「すべて展開」。フォルダ全体を展開してください。
3. 展開先の START.bat または GuidewireSimulator.exe をダブルクリック。

ZIPを開いた画面の中から直接実行せず、フォルダ全体を展開してください。
Windows 10/11の64bit版が対象です。exeは未署名です。
GitHubの「Code → Download ZIP」はソース版です。起動用ZIPはReleasesへ添付してください。

## 操作

- W / ↑：挿入、S / ↓：引き戻し。
- A / ←：左回転、D / →：右回転。キーを押している間だけ操作します。
- 設定 → Viewer設定 → 力学設定 → 剛性減衰βで変更。「適用」で保存。
- 減衰変更時は計算が再開始され、挿入状態がリセットされます。
- 下部の計算速度100%は実時間と同じ進行率。直近約1秒の速度を表示します。
- 設定とログは %LOCALAPPDATA%/GuidewireSimulator に保存されます。アプリの展開先へは書き込みません。

## この版の内容

2026-10-15作成のローカル版。NiTiの線形近似、半径0.18 mm、dt=0.002 s、
質量減衰α=0、剛性減衰β=0.02、摩擦係数0.15を初期値とします。
設定画面での減衰変更と、直近の実時間進行率表示を含みます。
人工的に生成したサンプル血管だけを同梱。患者形状・研究室データ・会話記録・個人PCの設定は含みません。
血管形状と全体表示の視野は公開用サンプルへ変更しています。計算核・材料値・時間刻みは元の試用版を保持しています。
研究室血管の操作結果がそのまま再現される配布物ではありません。

現時点で低減衰時のジャンプ・貫通や、多点接触時の計算遅延が観察されています。
実製品の反力・挙動は実測校正していません。動作確認と物理的な精度検証は別です。

## 新しい公開リポジトリへの置き方

- ソース版ZIPの中身を新しいリポジトリのルートへ置きます。
- Releasesを作り、起動用Windows ZIPとSHA256SUMS.txtを添付します。
- 配布ZIPには対応するソース・出典・第三者ライセンスを含めています。
- 既存の非公開リポジトリ、研究室データ、ローカルruntimeの実験ログはアップロードしないでください。

## ソースから実行・再ビルド

Python 3.9以降（Tk付き）で requirements.txt をインストール。
native/とupstream/はSOFA/BeamAdapter由来のコードを含み、出典はNOTICE.mdと各manifestに記録しています。
build_native.ps1 -FetchToolchain で計算核をビルド後、python portable_launcher.py で起動できます。
exe化には requirements-build.txt と build_portable.ps1 を使用します。
SOFA本体のインストールとGPUは不要です。

## 2026-10-15 起動修正

Error 13（portable_launcher.py 81行目）の対策として、ログ・設定・実行状態をユーザー用フォルダへ移しました。
複数プロセスで同じ起動ログを開かないようにしています。

## 問題があるとき

%LOCALAPPDATA%/GuidewireSimulator/startup-*.log と sessions/配下のsofa.logを確認してください。
元のフォルダがなくても動くよう、起動用exeにはPython・NumPy・Pillow・Tkを同梱しています。

## 出典・ライセンス

SOFA/BeamAdapter由来部分の著作権表示とLGPLを保持しています。
この公開スナップショットの自作接続部分もLGPL-2.1-or-laterとして配布します。
第三者ライブラリはそれぞれのライセンスに従います。LICENSE、NOTICE.md、THIRD_PARTY_NOTICES.md、licenses/を参照。
