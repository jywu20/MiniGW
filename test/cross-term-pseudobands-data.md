julia> transition_matrix_irreducible_1BZ_def(wfn, 120, 4000, 1, 1, 1693)
0.040302111809754015 - 0.03879153056138239im

julia> transition_matrix_irreducible_1BZ_def(wfn, 120, 3999, 1, 1, 1693)
-0.06427947218465854 + 0.048474126547253873im

julia> transition_matrix_irreducible_1BZ_def(wfn, 120, 3998, 1, 1, 1693)
0.024957860680757227 - 0.0019894256816043155im

julia> transition_matrix_irreducible_1BZ_def(wfn, 120, 3997, 1, 1, 1693)
0.012546497714126536 - 0.004954173968399836im

julia> transition_matrix_irreducible_1BZ_def(wfn, 119, 4000, 1, 1, 1693)
0.053834857933451885 - 0.05981516515280634im

julia> transition_matrix_irreducible_1BZ_def(wfn, 119, 3999, 1, 1, 1693)
0.04228172570390553 - 0.03683901509920638im

julia> transition_matrix_irreducible_1BZ_def(wfn, 119, 3998, 1, 1, 1693)
-0.013281237165684238 + 0.00047220042651992315im

julia> transition_matrix_irreducible_1BZ_def(wfn, 119, 3997, 1, 1, 1693)
0.023455640383064972 - 0.007903048142666633im


M = [0.040302111809754015 - 0.03879153056138239im  -0.06427947218465854 + 0.048474126547253873im  0.024957860680757227 - 0.0019894256816043155im  0.012546497714126536 - 0.004954173968399836im
0.053834857933451885 - 0.05981516515280634im  0.04228172570390553 - 0.03683901509920638im  -0.013281237165684238 + 0.00047220042651992315im  0.023455640383064972 - 0.007903048142666633im]

2×4 Matrix{ComplexF64}:
 0.0403021-0.0387915im  -0.0642795+0.0484741im   0.0249579-0.00198943im  0.0125465-0.00495417im
 0.0538349-0.0598152im   0.0422817-0.036839im   -0.0132812+0.0004722im   0.0234556-0.00790305im


With k_idx = 2, q_idx = 2:

4×4 Matrix{ComplexF64}:
    0.0126313+4.84201e-5im    -0.0038378-0.0112728im    0.00959319+0.00837052im   -0.0115179-0.0100703im
   -0.0108141+0.00483306im   -0.00872765-0.00934351im    0.0153059+0.000555291im   0.0127742+0.000654051im
 -0.000472784+0.0021597im    -0.00656669+0.0104698im    0.00237742+0.00182027im   0.00541022+0.00121971im
  -0.00132824-0.012458im    -0.000533717+0.00218401im  -0.00536106-0.00104829im    0.0028955-0.000623324im
With no pseudobands : 0.0017974692229075096 + 0.0im
With pseudobands    : 0.0007055308133726507 + 3.3881317890172014e-20im

```julia
n_range = 115:120
n′_range = 1997:2000
G_idx = 1693
k_idx = 1
q_idx = 1
```

```
6×4 Matrix{ComplexF64}:
 -0.00012628+0.000100455im    0.00218344+0.000545692im  -0.00249633-0.000547771im   -0.00649905-0.0101251im
  0.00202288+0.000895542im   -2.99902e-6+0.000211186im    -0.011026-0.00480755im    0.000940664+0.0023859im
 -0.00442077+0.00327955im    -0.00106098+0.00077932im    0.00117315-0.00103874im    -0.00129031-0.000601032im
  -0.0012829-0.000344867im     0.0053143+0.00153636im    -0.0014534+5.77312e-5im   -0.000632387-0.00133127im
  0.00022801-0.000103829im  -0.000792003-0.000109542im    0.0154426-0.00146874im     0.00646733+0.00596351im
 -0.00072251+0.000239234im  -0.000177248-6.55646e-5im    0.00871233-0.00129562im     -0.0119402-0.00989424im
With no pseudobands : 0.0010231349164562034 + 0.0im
With pseudobands    : 0.0008688225843863581 + 0.0im
```

```julia
n_range = 119:120
n′_range = 3999:4000
G_idx = 1527 
k_idx = 2
q_idx = 2
```

2×2 Matrix{ComplexF64}:
 0.000128246-0.000192557im  -0.00179157-0.00369385im
  0.00295771-0.00261697im    6.53654e-5+0.000124894im
With no pseudobands : 3.2524227004877903e-5 + 0.0im
With pseudobands    : 3.322023451503344e-5 + 0.0im

```julia
let cross = 0.0, direct = 0.0
    n′_1 = 2
    for n in eachindex(n_range)
        for n′_2 in eachindex(n′_range)
            if n′_2 == n′_1 
                continue
            end
            cross += M_nn′[n, n′_1]' * M_nn′[n, n′_2]
        end
    end

    for n in eachindex(n_range)
        direct += M_nn′[n, n′_1]' * M_nn′[n, n′_1]
    end
    
    println("Cross term sum : $cross")
    println("Direct term    : $direct")
end
```

```julia
let χ_pb = 0.0 
    direct_terms = []
    cross_terms = []
    for n in eachindex(n_range)
        for n′_1 in eachindex(n′_range)  
            for n′_2 in eachindex(n′_range)
                term = M_nn′[n, n′_1]' * M_nn′[n, n′_2]
                χ_pb += term
                if n′_1 == n′_2 
                    push!(direct_terms, term)
                else
                    push!(cross_terms, term)
                end
            end
        end
    end
    
    println("With pseudobands    : $χ_pb")
    println(direct_terms)
    println(cross_terms)
end
```

```julia
let χ_pb = 0.0,
    χ = 0.0

    M_nn′ = rand(ComplexF64, length(n_range), length(n′_range))

    for n in eachindex(n_range)
        for n′ in eachindex(n′_range)
            χ += M_nn′[n, n′]' * M_nn′[n, n′]
        end
    end

    println("With no pseudobands : $χ")

    for n in eachindex(n_range)
        for n′_1 in eachindex(n′_range)  
            for n′_2 in eachindex(n′_range)
                χ_pb += M_nn′[n, n′_1]' * M_nn′[n, n′_2]
            end
        end
    end
    
    println("With pseudobands    : $χ_pb")
end
```

```
julia> let χ_pb = 0.0,
           χ = 0.0
       
           M_nn′ = rand(ComplexF64, length(n_range), length(n′_range))
       
           for n in eachindex(n_range)
               for n′ in eachindex(n′_range)
                   χ += M_nn′[n, n′]' * M_nn′[n, n′]
               end
           end
       
           println("With no pseudobands : $χ")
       
           for n in eachindex(n_range)
               for n′_1 in eachindex(n′_range)  
                   for n′_2 in eachindex(n′_range)
                       χ_pb += M_nn′[n, n′_1]' * M_nn′[n, n′_2]
                   end
               end
           end
           
           println("With pseudobands    : $χ_pb")
       end
With no pseudobands : 19.174201450118176 + 0.0im
With pseudobands    : 64.02872715594225 - 5.551115123125783e-16im
```

```julia
n_range = 115:118
n′_range = 3999:4000
G_idx = 1527 
k_idx = 2
q_idx = 2
```

```
4×2 Matrix{ComplexF64}:
 0.000458199-0.0100222im   -0.000377924-0.0107472im
 -0.00940528+0.00520504im    0.00835395-0.00551256im
  0.00610207+0.00624864im    0.00280169-0.00214768im
 -0.00059551-0.00349095im    0.00869543-0.000523001im
With no pseudobands : 0.0006091977445264163 + 0.0im
With pseudobands    : 0.0006103923226445423 + 0.0im
```

```julia
n_range = 115:118
n′_range = 3999:4000
G_idx = 1680 
k_idx = 2
q_idx = 2
```

```
4×2 Matrix{ComplexF64}:
  -0.0161703-0.0111062im   -0.00287007+0.00929894im
  0.00650508-0.0072097im     0.0177866+0.00832203im
  -0.0117198-0.00720211im   0.00222948+0.00443065im
 -0.00466746+0.00168165im   -0.0134971-0.002977im
With no pseudobands : 0.0013889245629970092 + 0.0im
With pseudobands    : 0.0013865056874492476 + 0.0im
```

```julia
n_range = 115:118
n′_range = 139:142
G_idx = 1680 
k_idx = 2
q_idx = 3
```

```
4×4 Matrix{ComplexF64}:
  0.00217471+0.00243517im   0.00513529-0.000434053im   -0.00599057-0.00216243im   -0.0057718+0.00466522im
  0.00285483+0.00424245im  -0.00326389-0.000180049im   -0.00208605+0.00711832im  -0.00424974-0.00471251im
 -0.00533387-0.00385054im   0.00268904+0.00556194im   -0.000554593+0.00192583im    0.0104851+0.0130451im
  0.00544663+0.00295955im   0.00360202+0.00553416im    -0.00900153-0.0141044im    0.00203554+0.00010968im
With no pseudobands : 0.00099669905099698 + 0.0im
With pseudobands    : 0.0004934200737054574 - 1.0164395367051604e-20im
```

It seems we finally find a counterexample of the validity of pseudobands; 
but no worry, it might be possible that there is a G vector 
with which the transition matrix is much larger than 0.001, 
and for that G vector pseudobands is equivalent to the correct answer... 
but not quite:

```julia
n_range = 115:118
n′_range = 139:142
G_idx = 1 
k_idx = 2
q_idx = 3
```

```
4×4 Matrix{ComplexF64}:
 0.00254309+0.00180421im   0.00291685+0.00379275im   -0.00999053-0.0132357im    0.00228282+0.00272208im
 0.00485467-7.61584e-5im  -0.00300553-0.000871996im  -0.00343318-0.00109227im   -0.0160997-0.00423687im
  0.0203975+0.0160654im    -0.0111627-0.0189036im      0.0042936+0.007791im    -0.00924641-0.0130618im
 -0.0182858-0.0118932im    -0.0152599-0.021199im      0.00946957+0.0130114im    0.00591825+0.00645798im
With no pseudobands : 0.003628863277375637 + 0.0im
With pseudobands    : 0.000980729313415162 + 1.3552527156068805e-20im
```

Here `G_idx` is intentionally chosen to be `1`,
because the two sets of wave functions seem to peak in roughly the same region.

```julia
n_range = 115:118
n′_range = 997:1000
G_idx = 1 
k_idx = 2
q_idx = 3
```

```
4×4 Matrix{ComplexF64}:
  7.26491e-6+0.000126234im  -2.23178e-5-0.000180048im   1.88434e-5+5.07272e-5im    3.38123e-5-2.81549e-6im
 0.000142328-0.000112991im   9.28571e-5-8.49582e-5im    1.91214e-6-3.40557e-5im    5.16516e-5+1.67367e-5im
  5.01413e-5-6.739e-5im      3.40169e-5+9.30511e-5im    4.24777e-5-0.000109287im   2.85139e-5+0.00010086im
 -9.49959e-5-2.73914e-5im   -6.20084e-5+5.49566e-5im   -5.86209e-5-8.69782e-5im   -0.00011063-3.66323e-5im
With no pseudobands : 1.8878597092727466e-7 + 0.0im
With pseudobands    : 2.712008076147698e-7 + 1.6543612251060553e-24im
```

Of course we should try to find, again, what's the "center $\vb{G}$" 
of the high-energy band.

```julia
n_range = 115:118
n′_range = 997:1000
G_idx = 810 
k_idx = 2
q_idx = 3
```

```
4×4 Matrix{ComplexF64}:
  0.00370246+0.00579609im  -0.00294989+0.00127736im   -0.0117116+0.010579im    0.0106464+0.00253356im
  0.00134004+0.00299146im   0.00656549-0.00168654im  -0.00286666-0.0105854im     0.01022-0.0120506im
  0.00622591-0.00895016im   -0.0122106+0.00450369im   -0.0288596+0.00142121im  0.0114031+0.00900794im
 -0.00339154+0.0124587im   -0.00830151+0.00702801im   0.00021863-0.014594im    0.0193925-0.0214655im
With no pseudobands : 0.0035223225332717788 + 0.0im
With pseudobands    : 0.0020180483721636837 + 0.0im
```

```julia
n_range = 110:114
n′_range = 997:1000
G_idx = 810 
k_idx = 2
q_idx = 3
```

```
5×4 Matrix{ComplexF64}:
  -0.0169516+0.00366348im    -0.0133306-0.00180948im    0.0135665-0.0132397im   -0.00235473-0.0244131im
 0.000904676+0.0114552im     0.00791486+0.0124615im    -0.0198746+0.0116769im     -0.026375-0.0116478im
 -0.00901365-0.0117104im      0.0105021+0.00485171im   -0.0182311+0.0222034im     0.0230536+0.00151106im
  0.00425345+0.00194565im    0.00699246-0.00138914im  -0.00197209-0.00873737im   -0.0316478+0.00057503im
  0.00688638-0.000901689im  -0.00364237+0.00323949im   -0.0283467-0.014157im     0.00561538-0.00702883im
With no pseudobands : 0.007178337402167596 + 0.0im
With pseudobands    : 0.005237038856229398 + 5.421010862427522e-20im
```

```julia
n_range = [109]
n′_range = 997:1000
G_idx = 810 
k_idx = 2
q_idx = 3
```

```
1×4 Matrix{ComplexF64}:
 -0.00796404+0.0107267im  0.0144087-0.00942721im  -0.00164074+0.0245173im  -0.0139272-0.0128109im
With no pseudobands : 0.0014368451514649976 + 0.0im
With pseudobands    : 0.00025238584219276896 - 5.421010862427522e-20im
```

```julia
n_range = 117:120
n′_range = 3000:3002
G_idx = 810 
k_idx = 2
q_idx = 3
```

```
4×3 Matrix{ComplexF64}:
    0.0123608-0.00177468im  -0.00198338-0.00123553im  0.000178783-0.00357986im
 -0.000614521+0.00432229im  0.000274486+0.00365312im   0.00192742-0.0012396im
  -0.00270002+0.00890874im  -0.00129603+5.19089e-5im  0.000172303-0.00423892im
  -0.00639433+0.00777336im  -0.00401727+0.00130321im  0.000286773+0.00123452im
With no pseudobands : 0.00043907055742565935 + 0.0im
With pseudobands    : 0.00044849993234468484 + 0.0im
```

```julia
n_range = 117:120
n′_range = 1999:2002
G_idx = 810 
k_idx = 2
q_idx = 3
```

```
4×4 Matrix{ComplexF64}:
 -0.00910969-0.0227463im    0.017072+0.01687im       -0.019742+0.0119986im   -0.00267768-0.00399959im
   0.0238821-0.00307714im  0.0211559-0.0126041im   -0.00294934-0.00375596im   -0.0169857+0.0156771im
   0.0194643-0.0243892im   0.0501617-0.00877075im  -0.00207243+0.00877126im   0.00726709+0.0124667im
  -0.0281934+0.0424099im   0.0301961-0.00813274im  -0.00764728+0.0121265im   -0.00181993-0.00876277im
With no pseudobands : 0.011189987027120523 + 0.0im
With pseudobands    : 0.008070646302970757 + 1.0842021724855044e-19im
```