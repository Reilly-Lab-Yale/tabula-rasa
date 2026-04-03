# Paper Plan

## NB vs ZINB model selection results

Restricting to the condition each dataset would actually use in practice
(obs for datasets with a transfection reporter, CM for seelig which has none):

| Dataset   | Condition       | LRT Bonferroni sig | AIC prefers ZINB | mu shift NB->ZINB |
|-----------|-----------------|--------------------|------------------|--------------------|
| Seelig    | CM (mandatory)  | 0%                 | 0%               | 0.22%              |
| Cohen     | obs             | 54%                | 82.5%            | 2.87%              |
| Shendure  | obs             | 6.7%               | 25%              | 12.3%              |

### Conclusions

- **Seelig (CM)**: NB is the right model. Zero models significant, zero AIC
  wins, mu parameters identical to 4 decimal places.

- **Cohen (obs)**: ZINB is warranted. Majority significant under Bonferroni,
  mu stable (~3% shift).

- **Shendure (obs)**: NB is probably sufficient. Only 7% survive Bonferroni,
  BIC agrees (3% prefer ZINB), 12% mu shift.

### Collision rates

| Dataset   | Estimated collision % |
|-----------|-----------------------|
| Cohen     | 0.117%                |
| Shendure  | 0.053%                |
| Seelig    | 7.630%                |

---

## Figures

(to be filled in)

<!--
Example entry format:

### Figure N: [title]
- **Data source**: [which ortho/TSV/results file]
- **Plot type**: [violin, bar, heatmap, etc.]
- **Key message**: [what the reader should take away]
- **Script**: [path to notebook/script that generates it]
- **Status**: [planned / drafted / final]
-->

---

## Supplementary

(to be filled in)
