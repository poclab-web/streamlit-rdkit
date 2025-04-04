from rdkit import Chem
from rdkit.Chem import AllChem
from ase import Atoms
from ase.calculators.emt import EMT
from ase.vibrations import Vibrations
import matplotlib.pyplot as plt

# 1. RDKitで分子構造を生成・最適化
def generate_structure(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)  # 水素を付加
    AllChem.EmbedMolecule(mol)  # 3D座標生成
    AllChem.MMFFOptimizeMolecule(mol)  # 最適化
    return mol

# 2. RDKit構造をASEのAtomsオブジェクトに変換
def rdkit_to_ase(mol):
    conf = mol.GetConformer()
    positions = []
    symbols = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        positions.append([pos.x, pos.y, pos.z])
        symbols.append(atom.GetSymbol())
    return Atoms(symbols=symbols, positions=positions)

# 3. ASEで振動数計算を実行
def calculate_vibrations(atoms):
    atoms.set_calculator(EMT())  # 力場の設定（ここではEMT）
    vib = Vibrations(atoms)
    vib.run()
    frequencies = vib.get_frequencies()  # 振動数を取得
    vib.clean()  # 一時ファイルを削除
    return frequencies

# 4. 振動モードをプロット（波数が大きい方を左に）
def plot_vibrations(frequencies):
    frequencies = sorted(frequencies, reverse=True)  # 大きい順にソート
    plt.figure(figsize=(8, 6))
    plt.stem(frequencies, [1] * len(frequencies), basefmt=" ", use_line_collection=True)
    plt.gca().invert_xaxis()  # X軸を反転
    plt.xlabel("Wavenumber (cm⁻¹)")
    plt.ylabel("Relative Intensity (a.u.)")
    plt.title("IR Spectrum (Frequencies Only)")
    plt.grid()
    plt.show()

# 実行
smiles = "CCO"  # エタノール
mol = generate_structure(smiles)
atoms = rdkit_to_ase(mol)
frequencies = calculate_vibrations(atoms)
plot_vibrations(frequencies)
