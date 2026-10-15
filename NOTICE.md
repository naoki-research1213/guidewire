# 引き継いだソース

`upstream/sofa/`、`upstream/include/sofa/` とそれらから抽出したC++関数はSOFAの原ソース・著作権表示を保持しています。原ソースが宣言するLGPLの全文は `upstream/sofa/LICENSE-LGPL.md` にあります。

BeamAdapterの原ソース・抽出関数は原著作権表示を保持しています。原LICENSEとAuthorsは `upstream/BeamAdapter/LICENSE.md`、`upstream/BeamAdapter/Authors.md` にあります。出典のコミットは `docs/source_extraction.md` と各manifestに記録しています。

Viewer・診断処理は同じユーザー所有の現行リポジトリから取り込んでいます。取り込み元のコミット・ファイルハッシュと接続変更は `viewer/reference_manifest.json` と `tools/import_reference_ui.py` にあります。

原ソースの関数本体は保持し、Data・リンク・行列格納・所有権・起動の接続をこのリポジトリで実装しています。配布済みSOFA／BeamAdapterのDLLを本シミュレーターの計算核として再配布する構成ではありません。ビルドしたDLLとローカルのコンパイラ、患者メッシュ、生の比較状態はGit管理外です。
