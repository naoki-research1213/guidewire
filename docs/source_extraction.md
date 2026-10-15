# 現行コードの抽出・接続記録

現行の単一器具・固定STL血管の操作版を接続済み。原関数と接続部分を区別し、対象外経路・版の未確定事項は以下に残す。[起動とコードの読む順序](native_simulator.md)。

## 参照版

- 現行シーン: `guidewire---SOFA` の `9c9a013fcce03384620604381cd75c24aca1f531`。既存環境は変更しない。
- SOFA: 配布物 `v26.06.00/git-info.txt` に記載された `7c18e95d5c5f2839079892c69e7d89a313c79603`、公式 tag `v26.06.00` と一致。
- BeamAdapter: 公式 tag `release-v26.06`、`cac4005bd7c9f266c5c775ab88cebbd7c7fa83c7`。配布 DLL との厳密なビルド由来は未確認。版の同一性を確定した扱いにはしない。

## 抽出原則

`upstream/` に原文とライセンスを保存する。`tools/regenerate_native.py` が各extractスクリプトを呼び、保存済みの原ソースから関数本体をコピーし、クラス名・SOFA依存型の宣言部分を変更する。関数本体の変更有無とハッシュは `native/generated/*manifest.json` に記録する。SOFAの Vec/Mat/Quat/Transform/SpatialVector も原文を使う。

`native/adapters.h` の Data、アクセサ、状態格納、リンクはSOFAオブジェクトの接続を置き換える。数値処理との差を隠さず、各アダプタも検証対象とする。現行シーンで使わない機能は、無言で別処理へ置換しない。

## 入出力・数式・順序と移植先

| 処理 | 原実装 | 入力 → 出力・計算 | 移植先・未解決点 |
|---|---|---|---|
| サンプリング | IRC::computeInstrumentsCurvAbs / sortCurvAbs | xtip、区間長、nbBeams → 展開部材料座標。密度=nbBeams/区間長、先端から逆順の等間隔位置を生成、threshold=.01で重複除去 | beam_functions.inc、本体無変更。単一器具・剛体区間なしのシーン値をDataへ接続 |
| 挿入時の再配置 | IRC::totalLengthIsChanging / applyInterventionalRadiologyController | 新旧先端長差 ΔL → 補間対象座標 x−ΔL → 旧節点姿勢から新節点姿勢。要素基準長=x1−x0、器具内座標=x−(xtip−1000) | 同上。速度は原文が更新する範囲のみ。独自補間を加えない |
| 回転 | IRC::applyInterventionalRadiologyController / BaseBeamInterpolation::addBeam | rotationInstrument → 各要素のDOF-to-local変換の局所X軸Quat | addBeamアダプタは原処理対応。先端姿勢の直接指定は入っていない |
| 入口 | IRC::fixFirstNodesWithUntil / RestShapeSpringsForceField | 未使用節点の位置・静止位置・速度を固定。indexFirstNode=firstSimulatedNode−1、当該節点の並進・回転ばね=1e8 | IRC本体抽出済。原Rigidばねの力・線形化行列をkernel.hへ接続済み |
| 共回転・内力 | AdaptiveBeamForceFieldAndMass::addForce | 現姿勢 → スプライン中点フレーム、Adjoint、現形状と材料静止形状のSpatialVector差 u。f=−AᵀKu + 重力 | 原文抽出済。デフォルト reinforceLength=false、useShearStressComputation=true を保持 |
| 要素剛性・質量 | 同::computeStiffness / computeMass | E,ν,ρ,A,I,J,L → 各6×6ブロック。Asy=Asz=0、radius=.9 | 原文抽出済。単位・半径・物性を勝手に換算しない |
| 質量乗算 | 同::applyMassLarge / addMDx | 原文は入力datadxを使用せずm_gravityを使う | 原文を保持。疑わしい処理を修正して同一と称さない。この挙動も保持して実行DLLとの全体比較を実施 |
| 行列組立 | 同::addMBKToMatrix | M,K、Adjoint → AᵀMA、AᵀKA を節点ブロックに加算 | 元関数を抽出し、BTD行列格納へ接続済み |
| 陰的積分 | EulerImplicitSolver::solve | b=h[f−αMv+(h+β)K_force v]、A=(1+hα)M−h(h+β)K_force。Δv=A⁻¹b、v_free=v+Δv、q_free=q⊕hv_free | kernel.hで原ソースの演算順序に沿って接続。RigidCoordの指数写像は原文抽出 |
| 直接求解 | BTDLinearSolver::invert / solve / computeMinvBlock | 6×6三重対角ブロック → α⁻¹、λ、前進・後退代入、必要逆行列ブロック | btd_functions.inc、本体無変更。BTDMatrix::Blockも原文。subpartSolve=falseの初期値を保持。拘束用partial_solveと疎JMinvJtも抽出・接続済み |
| 接触用点 | MultiAdaptiveBeamMapping / IRC::interventionalRadiologyCollisionControls | 材料座標 → 元のスプラインによる点、J。衝突サンプル間隔 straight=19.6、spire=2 | IRCとMappingを接続済み。自由速度をJで写像し、自由接触点は現在位置+h×写像速度で更新 |
| 接触検出 | LocalMinDistance / BVHNarrowPhase | 三角形–点、辺–点、辺–辺。距離だけでなく局所最小フィルタと有効点判定。alarm=2、contact=.1、angleCone=.2 | 原関数・STL・トポロジー・BVH探索を抽出し接続済み。最近傍面やペナルティへの置換なし |
| 拘束・摩擦 | FrictionContact / BaseContactLagrangianConstraint / LCPConstraintSolver::nlcp_gaussseidel_unbuilt | 自由運動の拘束違反、3行接触基底、J → 法線相補・Coulomb摩擦、反復補正。μ=.15、tol=1e−4、maxIt=1000、build_lcp=false | LocalBlock33・unbuilt反復・補正・BTD・履歴を原文から抽出し、scene.hへ接続済み |
| 状態補正 | LinearSolverConstraintCorrection | 拘束力 → A⁻¹Jᵀλ、速度係数1、位置係数h（EulerImplicitSolverの定義）。既存ログはλ/h | 原関数をscene.hへ接続済み。位置係数h、速度係数1を保持 |
| 画面・入力 | 現行 xray_viewer / xray_physics / wire_raster | キー → 同じ挿入・回転コマンド、表示中心線 → 同じX線画像と計測表示 | 固定リビジョンから取り込み、抽出DLLワーカーへ接続済み。元診断メソッド本体を保持 |

全体のイベント・マッピング・衝突・補正の順序、warm start、検出ID、トポロジー順序を接続した。原CollisionResponseのモデル対ポインタ順に合わせ、現行で観測した三角形/点→辺/点→辺/辺の順序を保持する。材料境界ではWireRestShapeのEPSILONと選択条件を原文から抽出し、簡略分岐を廃止した。

標準・研究室血管の各2400ステップ、200ステップごとの12サンプルで、全節点姿勢・速度・衝突点・表示中心線・拘束インパルス・接触数を比較した。最大姿勢差は約1.3e−10、接触数は全サンプル一致。[数値記録](validation_runs/2026-10-11.json)。これは上記の条件における数値比較であり、全条件の同一性や実物の精度を証明するものではない。

## 接続途中の記録（以下は当時の状態）

以下の未接続・未完成という記述は当時の状態を記録したもの。現在の対応は上表と数値記録を参照する。

## 2026-10-11 接続確認

接触なしの `Kernel::advanceWithoutContact` を、元シーンの固定血管を初期化前に遠方へ移動した比較実行と照合した。元リポジトリのファイルは変更していない。h=.005、400step、挿入軸+1、201stepまでは回転0、202step以降は回転軸+1。40step毎の10点における先端中心位置差の最大絶対値は2.1540103034567437e−11。最終xtip=20.00000011920944、先端=(-10.0626290539,4.45264409451,-12.3740319211)。この比較は接触・長時間・全節点・性能の保証ではない。

BTDの部分求解、疎JMinvJt、LinearSolverConstraintCorrectionのresetForUnbuiltResolution / addConstraintDisplacement / setConstraintDForce / getBlockDiagonalCompliance / applyContactForce、LocalBlock33のcompute / New_GS_State、LCPConstraintSolverのnlcp_gaussseidel_unbuiltを原文で抽出した。3行の合成拘束による内部確認で4反復・残差6.07e−5、法線非負とμ=.15の摩擦円錐を満たす。これは血管接触の同一性確認の代わりにはならない。

Kernelの陰的積分はSOFAのシーングラフ演算を固定データ配列への呼び出しに展開しており、EulerImplicitSolver::solve全文の無変更コピーではない。原ForceFieldのaddForce / addMDx / applyStiffnessLarge / addMBKToMatrix、原Spring::addForce、原BTD::invert / solveを同じ順序で呼ぶ。固定条件の射影、入口ばねの線形化行列とRayleigh項の接続も対応付けた。上記比較の範囲で数値一致を確認したが、未接続のマッピング幾何剛性・接触履歴を含めた差異はまだ検証が必要。

### 接触幾何・写像の抽出

LocalMinDistanceの三角形–点／辺–点／辺–辺と3種類のtestValidity、CubeCollisionModelの境界箱・安定ソートによる分割、BVHNarrowPhaseの探索、Point/Line/TriangleCollisionModelのcomputeBoundingTreeを関数本体無変更で抽出した。MeshSTLLoaderのbinary/ASCII読み込みとMeshTopologyの辺・隣接配列生成も原文を保持。geometry_probeで標準血管（1274頂点、2516三角形、3790辺）と研究室血管（9792頂点、19520三角形、29312辺）の読み込み・探索・接触検出が通った。局所的な合成点の確認であり、シーン全体の一致確認ではない。

AdaptiveBeamMappingのcomputeDistribution/apply/applyJTonPoint/拘束applyJT、BarycentricMapperMeshTopologyの重心座標計算・apply/拘束applyJT、BaseContactLagrangianConstraintの接触基底・拘束行・自由状態違反、BaseUnilateralContactResponseの重複除去・mapper接続・response生成も本体無変更で抽出した。Data、疎行列、モデル所有権、逐次task実行はアダプタで接続する。現在のシーンで固定のBVH深さ6を対象とし、実行中の深さ変更は未対応として明示的に停止する。現在シーンで設定されていないCollisionElementActiverを独自に追加しない。

接触の写像と全体ステップへの接続・履歴・表示への出力は引き続き作業中。まだ完成版として起動入口を切り替えていない。
