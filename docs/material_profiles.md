# 公開材料物性を反映した比較版

試用設定の数値補償とは別に、材料メーカーが公表するEと密度を直接入力する比較版を追加した。実製品のCHIKAI/Synchro等を実測して校正した設定ではない。

| 起動ファイル | 材料 | E [Pa] | ρ [kg/m³] |
|---|---|---|---|
| `launch_native_oi_steel316lvm.bat` | 316LVM鋼 | 2.00×10¹¹ | 8000 |
| `launch_native_oi_nitinol.bat` | NiTiオーステナイト、線形近似 | 7.50×10¹⁰ | 6500 |

鋼の値は20°Cの[Alleima 316LVM材料資料](https://www.alleima.com/en/technical-center/material-datasheets/bar-and-hollow-bar/bar/Alleima-316lvm/)。NiTiは[Jacques Allemann JAS-036、2026/03版](https://www.jacques-allemann.com/file/2478/Legierungen%20/%20Tech-Sheets/JAS%20-%20Tech%20Sheet%20/%20Werkstoffe/jas-036-tech-sheet-nitinol.pdf)のオーステナイトE=75–80 GPaの下限を選んだ。対象製品の鋼種・NiTiの相状態を特定したものではない。

## 設定の意味

- どちらも直径0.36 mmの均一な中実円断面。テーパーしたコア、コイル、被覆の複合構造を再現していない。
- 半径、表示直径、接触離隔、模擬する区間長、自然曲線、操作速度は先の実寸試用と同じ。CHIKAIの全長2 mへ延長していない。
- ν=0.3は仮定。μ=0.15、剛性減衰β=0.2 s、入口拘束ばねは数値設定であり、対象製品の測定値ではない。
- 実物E・ρへ合わせるため、以前のEI/GJ/線質量を保存する補償は使わない。中実仮定のEIは鋼1.64896×10⁻⁴、NiTi6.18360×10⁻⁵ N·m²。線密度は鋼8.14301×10⁻⁴、NiTi6.61619×10⁻⁴ kg/m。
- NiTiは一定Eの線形材料として使う。相変態・超弾性の応力平坦部・履歴は原モデルに追加していない。

操作して硬く感じたり反力が大きくなったりしても、実製品が同じ挙動になると判断しない。部位別の完成品EI/GJ・線密度を取得し、等価モデルを同定する作業が残る。

## SI値を編集する

設定は`config/guidewire_material_steel316lvm.json`と`config/guidewire_material_nitinol.json`。`material_si`のEはPa、密度はkg/m³、重力はm/s²。`boundary_si`の入口並進ばねはN/m、回転ばねはN·m/rad。編集後にViewerを再起動する。再ビルドは不要。

長さ単位ℓ₀=1/6000 m、質量単位m₀=1 kg、時間単位1 sを明示し、既存計算座標へ次の換算を行う。

\[
E_{code}=E_{SI}\ell_0/m_0,\quad
\rho_{code}=\rho_{SI}\ell_0^3/m_0,\quad
g_{code}=g_{SI}/\ell_0,
\]
\[
k_{code}=k_{SI}/m_0,\quad
k_{\theta,code}=k_{\theta,SI}/(m_0\ell_0^2),\quad
F_{SI}=F_{code}m_0\ell_0.
\]

重力は−9.81 m/s²、コードでは−58860。接触反力と並進支持反力は、従来の再構成式で求めた結果をNへ換算して表示・CSV保存する。これは単位の校正であり、モデルの実験的な精度保証ではない。ログの座標は既存の計算座標のまま保持し、JSON/CSVに長さ尺度を明記する。支持反力は入口ばねとその剛性減衰の反力であり、総挿入駆動力そのものではない。

原実装の`addMDx`は速度入力を使わず重力を読むため、この比較版だけ質量Rayleigh係数α=0にした。非標準の質量減衰項を係数で無効化し、原関数は変更していない。既定版・以前の実寸試用はα=0.1を保持する。

## 確認した範囲

各材料・各血管で6000ステップ、12秒分の挿入・回転・引戻しを実施。有限状態、単位Quat、接触発生、引戻しを確認した。鋼の初回確認は反復上限1000回で2ステップ未収束だったため、許容誤差1e−4を変えず反復上限を5000回に変更した。最終設定は4ケースとも未収束0。[確認記録](validation_runs/2026-10-11-material-profiles.json)。

質量の数値単位を1→0.01 kgとして1500ステップ比較し、姿勢差最大約6.0×10⁻¹¹、N換算した拘束乗数差最大約4.6×10⁻⁸を確認した。[単位換算の比較](validation_runs/2026-10-11-material-units.json)。

両材料とも研究室Viewerの操作コールバックで挿入・回転・引戻し、投影画像、N表示の出力を確認した。人がキーを操作して確認した結果や、実験との一致を証明する結果ではない。計算核のみの12秒分の実行は約5–8.5秒であり、Viewer・ログを含む実時間性とは分けて扱う。

再確認は`tools/check_material_profiles.py`、単位の比較は`tools/check_material_units.py`。起動設定の作成処理は`tools/create_material_profiles.py`で、ユーザーが編集した設定を再生成すると上書きするため、試す設定は別名で複製できる。
