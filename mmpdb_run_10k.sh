#!/usr/bin/env bash
mmpdb fragment compound_smiles_10k.csv -o compound_smiles_10k.fragments
mmpdb index compound_smiles_10k.fragments -o compound_smiles_10k.mmpdb