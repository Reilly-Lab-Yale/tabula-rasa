# Bug: Bounds.from_ortho() fails on seelig_ortho_20260320

## Symptom

```
TypeError: float() argument must be a string or a real number, not 'Series'
```

Traceback points to `core.py` in `Bounds.from_ortho()`:

```python
for key in getattr(inp, var).nb:
    current = getattr(inp, var).nb[key].result()
    maxes.append(float(current.max()))   # <-- crashes here
```

`current.max()` returns a `Series` (one max per column) when `current` is a
DataFrame, but `float()` expects a scalar.

## Root cause (suspected)

The NB parameter storage format changed at some point during ortho refactoring.
Earlier orthos (shendure_ortho_20260306, cohen_ortho) store `nb[key]` as a
1-D structure whose `.result()` returns a scalar or Series that reduces to a
scalar under `.max()`. The seelig_ortho_20260320 stores it as a 2-D DataFrame,
so `.max()` returns a column-wise Series instead.

## Impact

- `Bounds.from_ortho()` cannot be run on seelig_ortho_20260320.
- Existing shendure and cohen bounds presets are unaffected (already saved).
- No seelig_bounds.tgz preset exists yet.

## Workaround used

Collision rate was computed directly from the raw TSV
(`seelig_collision_rate.py`) without going through `from_ortho`. The value
7.630123998189407% was recorded in `estimated_percent_conflict.ipynb`.

## Fix (deferred)

In `Bounds.from_ortho()` around the min/max loop, replace:

```python
maxes.append(float(current.max()))
mins.append(float(current.min()))
```

with something robust to both 1-D and 2-D results, e.g.:

```python
maxes.append(float(np.nanmax(current.values if hasattr(current, 'values') else current)))
mins.append(float(np.nanmin(current.values if hasattr(current, 'values') else current)))
```

but first confirm what the new parameter storage format actually looks like
before patching, to avoid masking a deeper structural mismatch.
