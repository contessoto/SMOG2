# Quickstart

This quickstart uses 2CI2 because it is small and works well for both all-atom
and C-alpha examples.

## 1. Download A PDB

```bash
curl -fsSLO https://files.rcsb.org/download/2CI2.pdb
```

If you are working from this source checkout without network access, use:

```bash
cp SMOG-CHECK/share/PDB.files/2ci2_v2.pdb 2CI2.pdb
```

## 2. Generate An All-Atom Model

```bash
smog3 -i 2CI2.pdb -AA \
  -o model.top \
  -g model.gro \
  -n model.ndx \
  -s model.contacts
```

Expected outputs:

- `model.top`
- `model.gro`
- `model.ndx`
- `model.contacts`

## 3. Generate A C-Alpha Model

```bash
smog3 -i 2CI2.pdb -CA \
  -o model_ca.top \
  -g model_ca.gro \
  -n model_ca.ndx \
  -s model_ca.contacts
```

## 4. Generate OpenSMOG XML

The OpenSMOG XML flag spelling is:

```bash
-OpenSMOG -OpenSMOGxml model.xml
```

Example:

```bash
smog3 -i 2CI2.pdb -AA -OpenSMOG -OpenSMOGxml model.xml \
  -o model_os.top \
  -g model_os.gro \
  -n model_os.ndx \
  -s model_os.contacts
```

## 5. Inspect Outputs

```bash
head -20 model.top
head -5 model.gro
wc -l model.contacts
```

The generated model files are local artifacts and should not be committed.
