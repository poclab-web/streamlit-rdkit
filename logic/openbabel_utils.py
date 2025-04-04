from openbabel import pybel

# SMILESから3次元構造を生成しエネルギー最適化する関数
def smiles_to_3d_with_make3D(smiles, forcefield='mmff94', steps=1000):
    # Pybel分子オブジェクトを作成
    mol = pybel.readstring("smi", smiles)

    # 3次元構造の生成とエネルギー最適化
    mol.make3D(forcefield=forcefield, steps=steps)

    # 3次元構造をMOL形式で返す
    return mol.write("mol")

# SMILESから3次元構造を予測
smiles = "CCO"  # エタノール
mol_3d = smiles_to_3d_with_make3D(smiles)

print("生成された3次元構造（MOL形式）:")
print(mol_3d)
