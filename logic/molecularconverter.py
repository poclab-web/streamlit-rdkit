import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem, SDWriter
from rdkit.Geometry import Point3D


class MoleculeConverter:
    """
    分子情報をSMILES、XYZ、SDF間で相互変換し、Z-matrixを生成するクラス。
    """
    def __init__(self):
        self.mol = None
        self.conf = None

    def load_from_smiles(self, smiles):
        """SMILESから分子を生成し3D座標を付与する"""
        self.mol = Chem.MolFromSmiles(smiles)
        if self.mol:
            self.mol = Chem.AddHs(self.mol)  # 水素原子を追加
            AllChem.EmbedMolecule(self.mol, AllChem.ETKDG())  # 3D座標生成
            self.conf = self.mol.GetConformer()
        else:
            raise ValueError("無効なSMILESです。")

    def load_from_xyz(self, xyz_coordinates):
        """XYZ座標情報を元に分子オブジェクトを生成する"""
        mol = Chem.RWMol()
        for element, _, _, _ in xyz_coordinates:
            mol.AddAtom(Chem.Atom(element))
        self.mol = mol.GetMol()
        conf = Chem.Conformer(len(xyz_coordinates))
        for i, (_, x, y, z) in enumerate(xyz_coordinates):
            conf.SetAtomPosition(i, Point3D(x, y, z))
        self.mol.AddConformer(conf)
        self.conf = self.mol.GetConformer()

    def load_from_sdf(self, sdf_file):
        """SDFファイルから分子を読み取る"""
        suppl = Chem.SDMolSupplier(sdf_file)
        if not suppl or suppl[0] is None:
            raise ValueError("無効なSDFファイルです。")
        self.mol = suppl[0]
        self.conf = self.mol.GetConformer()

    def export_xyz(self):
        """分子のXYZ形式を出力"""
        if not self.conf:
            raise ValueError("分子データがありません。")
        output = []
        for atom in self.mol.GetAtoms():
            pos = self.conf.GetAtomPosition(atom.GetIdx())
            output.append(f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}")
        return "\n".join(output)

    def export_sdf(self, output_filename):
        """SDF形式のファイルに分子をエクスポート"""
        if not self.mol:
            raise ValueError("分子データがありません。")
        writer = SDWriter(output_filename)
        writer.write(self.mol)
        writer.close()

    def generate_zmatrix(self):
        """分子のZ-matrixを構築する"""
        if not self.conf:
            raise ValueError("分子データがありません。")
        zmatrix = []
        atoms = self.mol.GetAtoms()

        def distance(i, j):
            """2つの原子間の距離を計算"""
            return np.linalg.norm(
                np.array(self.conf.GetAtomPosition(i)) - np.array(self.conf.GetAtomPosition(j))
            )

        def angle(i, j, k):
            """3点間の結合角を計算"""
            p1, p2, p3 = [self.conf.GetAtomPosition(x) for x in (i, j, k)]
            v1 = np.array([p1.x - p2.x, p1.y - p2.y, p1.z - p2.z])
            v2 = np.array([p3.x - p2.x, p3.y - p2.y, p3.z - p2.z])
            cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
            return np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))

        def dihedral(i, j, k, l):
            """4点間の二面角を計算"""
            p1, p2, p3, p4 = [self.conf.GetAtomPosition(x) for x in (i, j, k, l)]
            b1 = np.array([p2.x - p1.x, p2.y - p1.y, p2.z - p1.z])
            b2 = np.array([p3.x - p2.x, p3.y - p2.y, p3.z - p2.z])
            b3 = np.array([p4.x - p3.x, p4.y - p3.y, p4.z - p3.z])
            n1 = np.cross(b1, b2)
            n2 = np.cross(b2, b3)
            m1 = np.cross(n1, n2)
            x = np.dot(n1, n2)
            y = np.dot(m1, b2) * np.linalg.norm(b2)
            return np.degrees(np.arctan2(y, x))

        for i, atom in enumerate(atoms):
            if i == 0:
                zmatrix.append((atom.GetSymbol(), None, None, None, None, None, None))
            elif i == 1:
                zmatrix.append((atom.GetSymbol(), 0, distance(0, 1), None, None, None, None))
            elif i == 2:
                zmatrix.append((atom.GetSymbol(), 1, distance(1, 2), 0, angle(0, 1, 2), None, None))
            else:
                zmatrix.append((
                    atom.GetSymbol(), 2, distance(2, i),
                    1, angle(1, 2, i),
                    0, dihedral(0, 1, 2, i)
                ))
        return zmatrix
