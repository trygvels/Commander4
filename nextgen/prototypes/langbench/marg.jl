using Random, Base.Threads
const npix = 12*512^2; const nb = 8
const f = [23,30,44,70,100,143,217,353.0]; const hk = 0.04799243
g = @. (hk*f/2.7255)^2*exp(hk*f/2.7255)/expm1(hk*f/2.7255)^2
mix(bs,bd) = hcat(g, (f./30).^bs, (f./353).^(bd+1) .* expm1(hk*353/20) ./ expm1.(hk.*f./20))
function marg(A, W, d)
    np = size(d, 2)
    parts = zeros(nthreads())
    @threads :static for t in 1:nthreads()
        s = 0.0
        for p in t:nthreads():np
            M00=M01=M02=M11=M12=M22=b0=b1=b2=0.0
            @inbounds @simd for k in 1:nb
                wk=W[k,p]; x0=A[k,1]; x1=A[k,2]; x2=A[k,3]; wd=wk*d[k,p]
                M00+=wk*x0*x0; M01+=wk*x0*x1; M02+=wk*x0*x2; M11+=wk*x1*x1; M12+=wk*x1*x2; M22+=wk*x2*x2
                b0+=wd*x0; b1+=wd*x1; b2+=wd*x2
            end
            c00=M11*M22-M12*M12; c01=M02*M12-M01*M22; c02=M01*M12-M02*M11
            c11=M00*M22-M02*M02; c12=M01*M02-M00*M12; c22=M00*M11-M01*M01
            det=M00*c00+M01*c01+M02*c02
            q=(b0*(c00*b0+c01*b1+c02*b2)+b1*(c01*b0+c11*b1+c12*b2)+b2*(c02*b0+c12*b1+c22*b2))/det
            s += 0.5q - 0.5log(det)
        end
        parts[t] = s
    end
    sum(parts)
end
rng = Xoshiro(1); w = 1 ./ [4,3,3,3.5,1.2,.8,1.2,4.0].^2
a = vcat(70randn(rng,1,npix), 10 .^ (2.5rand(rng,1,npix)), 10 .^ (0.5 .+ 2.5rand(rng,1,npix)))
d = mix(-3.05,1.58)*a .+ randn(rng,nb,npix)./sqrt.(w)
W = repeat(w, 1, npix)
println("check ", marg(mix(-3.0,1.6), W, d))
t = @elapsed for i in 0:4; marg(mix(-3.0+0.01i,1.6), W, d); end
println("julia ($(nthreads()) threads): ", round(t/5, digits=3), " s")
