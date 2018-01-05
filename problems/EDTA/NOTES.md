The definition 

```
    var moves : [0..n, 0..m] move;
```

in `edit_alignment` fails mysteriously at compilation if the enum
is defined within that procedure rather than at module scope.
