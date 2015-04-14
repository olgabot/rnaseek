# An assortment of helpful commands

## Check if MISO ran successfully

From the sample directory, run this command, replacing `M2nd_33` with the sample you care about. Broken down, this is what this command is doing:

1. `ls -lh miso/M2nd_33*/*/summary/*`: Searches for all `*.miso_summary` files
2. `grep -v 113`: Removes the ones of size 113 bytes, which is exactly the size of the MISO header (will break down if `113` is anywhere else in the filename)
3. `grep -v ' 0 '`: Removes files of size 0
4. `wc -l`: Counts the remaining lines. If this is equal to 8, then all of the MISO splice types have been calculated

```
ls -lh miso/M2nd_33*/*/summary/* | grep -v 113 | grep -v ' 0 ' | wc -l
```