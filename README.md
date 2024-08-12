# DiscreteDiff

[![Build Status](https://github.com/tkrhsmt/DiscreteDiff.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tkrhsmt/DiscreteDiff.jl/actions/workflows/CI.yml?query=branch%3Amain)

DiscreteDiffは，等分割格子における有限差分法を利用して，3次元配列の導関数を求めます．

## インストール手順

`DiscreteDiff.jl`をインストールするには，juliaのREPLモードでパッケージモードに変更します．

```julia
]
```

この後，次のインストールコマンドを利用します．

```julia
add https://github.com/tkrhsmt/DiscreteDiff.jl
```

## 使用できる差分スキーム

`DiscreteDiff.jl`では，以下に示された複数の差分スキームを使用できます．

- 1階中心差分
    - 2次精度
    - 4次精度
    - 6次精度
- 1階compact差分
    - 6次精度
- 中心補間
    - 2次精度
    - 4次精度
    - 6次精度
- compact補間
    - 6次精度
- 1階中心補間差分
    - 2次精度
    - 4次精度
    - 6次精度
- 1階compact補間差分
    - 6次精度
- 2階中心差分
    - 2次精度
    - 4次精度
    - 6次精度
- 2階compact差分
    - 6次精度

全てのスキームに対して，以下の境界条件を指定できます．

- ディリクレ境界(値を指定可能)
- ノイマン境界(0のみ)
- 周期境界

## 使用できる関数

上記の差分スキームを利用するには，以下の関数を利用します．

### 1階微分
```julia
d1x(u, dx, scheme,  boundary, dirichlet)
```

- `u` : 微分を行う3次元配列
- `dx` : 差分幅 (float型)
- `scheme` : スキームの選択 (int型)
    - `scheme = 1` : 2次精度中心差分
    - `scheme = 2` : 4次精度中心差分
    - `scheme = 3` : 6次精度中心差分
    - `scheme = 4` : 6次精度compact差分
- `boundary` : 境界条件の指定 (int型 1次元6成分配列)
    - `0` : 周期境界
    - `1` : ノイマン境界
    - `2` : ディリクレ境界
- `dirichlet` : ディリクレ条件の値 (float型 1次元6成分配列)

> - `boundary`における境界指定では，1次元6成分配列を入力する必要があります．
> このとき， `[<x方向境界1つ目>,<x方向境界2つ目>,<y方向境界1つ目>,<y方向境界2つ目>,<z方向境界1つ目>,<z方向境界2つ目> ]`となります．
> - `dirichlet`におけるディリクレ条件の指定は，1次元6成分配列を入力する必要があります．
> このとき，`[<x方向境界1つ目>,<x方向境界2つ目>,<y方向境界1つ目>,<y方向境界2つ目>,<z方向境界1つ目>,<z方向境界2つ目> ]`となります．
> - `dirichlet`における指定は，`boundary`における指定でディリクレ境界を選んだ場合のみに使用されます．
> - `scheme,　boundary,　dirichlet`は， デフォルト値が設定されています．デフォルトでは，2次中心差分，3方向周期境界，3方向0のディリクレ境界です．

同様にすることで，y方向の差分に`d1y`，z方向の差分に`d1z`関数が使用できます．

### 2階微分

2階微分でも，1階微分と同じ形で記述できます．

```julia
d2x(u, dx, scheme,  boundary, dirichlet)
```

それぞれの項目は，1階微分と同じです．
同様にして，y方向の差分に`d2y`，z方向の差分に`d2z`関数が使用できます．

### 中心補間

```julia
intx(u, scheme,  boundary, dirichlet, mode)
```

それぞれの項目は，1階微分と同じです．
`dx`が省略されていますが，補間場所は全て中心が使用されるためです．
`mode`は前進補間，および後退補間を選択できます．
- `1` : 前進補間
- `2` : 後退補間

同様にして，y方向に`inty`，z方向に`intz`関数が使用できます．

### 1階中心補間差分

```julia
int_d1x(u, dx, scheme,  boundary, dirichlet, mode)
```

それぞれの項目は，1階微分と同じです．
`mode`は前進補間，および後退補間を選択できます．
- `1` : 前進補間
- `2` : 後退補間

同様にして，y方向の差分に`int_d1y`，z方向の差分に`int_d1z`関数が使用できます．
