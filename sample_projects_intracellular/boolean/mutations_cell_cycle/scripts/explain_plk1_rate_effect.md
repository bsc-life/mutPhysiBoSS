# Why PLK1 Rate Mutations Don't Correlate with PLK1 Node State

## The PLK1 Node Logic

Looking at the Boolean network definition:

```bnd
Node Plk1 {
  logic = (!FoxM1 & !Cdc25A & CyclinB & Cdk1 & !Cdh1 & Plk1_H) | 
          (FoxM1 & !Cdc25A & CyclinB & Cdk1 & !Cdh1) | 
          ... [other conditions];
  rate_up = @logic ? $u_Plk1 : 0;
  rate_down = @logic ? 0 : $d_Plk1;
}
```

## How MaBoSS Rate Parameters Work

**Critical Understanding:**
- When `logic = TRUE`: `rate_up = $u_Plk1` and `rate_down = 0` → PLK1 turns **ON**
- When `logic = FALSE`: `rate_up = 0` and `rate_down = $d_Plk1` → PLK1 turns **OFF**

## Why There's No Correlation

1. **The Boolean logic determines the DIRECTION** (ON vs OFF)
   - If upstream nodes (CyclinB, Cdk1, FoxM1, etc.) are in the right state, the logic is TRUE
   - When logic is TRUE, PLK1 will turn ON **regardless of the down rate value**
   - The down rate is set to 0 when logic is TRUE, so it has no effect

2. **The rate parameter only affects the SPEED** of transitions
   - A higher `$d_Plk1` makes PLK1 turn OFF **faster** when logic is FALSE
   - But it doesn't change **when** PLK1 turns OFF (that's determined by the logic)

3. **The logic is the dominant factor**
   - During cell cycle phases where CyclinB, Cdk1, etc. are active, PLK1 logic is TRUE
   - PLK1 will be ON during these phases regardless of mutation status
   - The mutation only affects how quickly PLK1 turns OFF when the logic becomes FALSE

## What This Means

- **Increasing `$d_Plk1`** (down rate) makes PLK1 turn OFF faster when conditions are not met
- **But it doesn't prevent PLK1 from turning ON** when the Boolean logic conditions are satisfied
- The **steady-state probability** of PLK1 being ON is primarily determined by:
  - How often the Boolean logic is TRUE (cell cycle phase distribution)
  - The ratio of up_rate to down_rate when logic is TRUE vs FALSE
  
## Why We See Low PLK1 ON Percentage

The low PLK1 ON percentage (0.34% average) suggests that:
- Most cells are in phases where the PLK1 activation logic is FALSE
- The Boolean network naturally keeps PLK1 OFF most of the time
- Mutations to the down rate don't change this because the logic is the controlling factor

## Potential Solutions

If you want mutations to affect PLK1 state more directly, you could:

1. **Mutate the upstream nodes** that control PLK1 logic (e.g., CyclinB, Cdk1, FoxM1)
2. **Mutate the up rate** (`$u_Plk1`) instead of the down rate
3. **Mutate the Boolean logic itself** (more complex, requires network modification)
4. **Use a much larger effect size** to see subtle effects on transition timing

