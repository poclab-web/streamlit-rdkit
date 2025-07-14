from collections import Counter
import selfies as sf

def count_selfies_features(smiles):
    try:
        tokens = sf.split_selfies(sf.encoder(smiles))
    except:
        return None

    counts = Counter(tokens)

    def count_contains(substr):
        return sum(v for k, v in counts.items() if substr in k)

    # 特徴量辞書
    features = {
        # 原子種カウント（no.1〜15）
        "C_aliphatic": counts["[C]"],
        "C_aromatic": counts["[c]"],
        "N_aliphatic": counts["[N]"],
        "N_aromatic": counts["[n]"],
        "O_aliphatic": counts["[O]"],
        "O_aromatic": counts["[o]"],
        "S": counts["[S]"] + counts["[s]"],
        "P": counts["[P]"] + counts["[p]"],
        "F": counts["[F]"],
        "Cl": counts["[Cl]"],
        "Br": counts["[Br]"],
        "I": counts["[I]"],
        "B": counts["[B]"],

        # その他の原子（no.16）
        "special_atoms": sum(
            v for k, v in counts.items()
            if k.startswith("[") and k not in {
                "[C]", "[c]", "[N]", "[n]", "[O]", "[o]", "[S]", "[s]",
                "[P]", "[p]", "[F]", "[Cl]", "[Br]", "[I]", "[B]",
                "[H]", "[Ring1]", "[Ring2]", "[Ring3]", "[Ring4]", "[Ring5]",
                "[Ring6]", "[Ring7]", "[Ring8]", "[Ring9]",
                "[Branch1]", "[Branch2]"
            }
        ),

        # 単結合（no.17） → SELFIESでは省略されるため記録しない

        # 二重・三重結合（no.18, 19）
        "double_bonds": sum(v for k, v in counts.items() if "=" in k),
        "triple_bonds": sum(v for k, v in counts.items() if "#" in k),

        # 特殊記号（[ ]）はSELFIESには現れない（no.20）

        # 電荷（no.21, 22）
        "positive_charge": count_contains("+1"),
        "negative_charge": count_contains("-1"),

        # 水素（no.23）
        "explicit_H": counts["[H]"],

        # 分岐（no.24）
        "branches": counts["[Branch1]"] + counts["[Branch2]"],

        # 環（no.25〜34）
        "rings": sum(counts[t] for t in counts if "Ring" in t),
        "ring_types": len(set(t for t in tokens if "Ring" in t)),

        # 立体中心（補足：記号に"@"が含まれるトークン）
        "stereo_centers": count_contains("@"),
    }

    return features
