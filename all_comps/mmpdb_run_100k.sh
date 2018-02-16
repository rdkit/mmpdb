#!/usr/bin/env bash
mmpdb fragment compound_smiles_100k.csv -o compound_smiles_100k.fragments
mmpdb index compound_smiles_100k.fragments -o compound_smiles_100k.mmpdb
