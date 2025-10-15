# Examples

## Prerequisities

Nauty binaries must be available. Note that nauty libraries are required to run the apps in the first order.

> [!NOTE]
> In the examples we will fix the number of vertices/dimension of real Bott manifolds to 8. This should give quite straightforward picture in an efficient way of how the apps work.

> [!WARNING]
> We assume that applications are run from the project's folder. In the folder we have created a `data` subfolder.

## Generating real Bott manifolds

```bash
$ geng 8 | directg -a | ./orbitg -o data/8-bott.d6
```

Alternatively

```bash
$ geng 8 | directg -a | ./minimalf -o data/8-bott.d6
```

## Orientable real Bott manifolds

```bash
$ cat data/8-bott.d6 | ./orientedf > data/8-bott-orientable.d6
```

## Spinc real Bott manifolds

```bash
$ cat data/8-bott-orientable.d6 | ./spincf > data/8-bott-spinc.d6
```

Alternatively

```bash
$ ./mats -d 8 | ./canonicalg | ./uniqueg > data/8-bott-spinc.d6
```

Alternatively

```bash
$ ./backtrack -d 8 | ./canonicalg | ./uniqueg > data/8-bott-spinc.d6
```

## Spin real Bott manifolds

```bash
$ cat data/8-bott-spinc.d6 | ./spinf > data/8-bott-spin.d6
```