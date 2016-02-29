function recvec(v; compute=BigInt, output=Int64)
  n = length(v)
  l = Array(Rational{compute}, n)
  l[:] = v[:]
  sl = n
  z = div(n, 2)
  q1 = Array(Rational{compute}, n)
  q2 = Array(Rational{compute}, n)
  q1[1] = 0
  q2[1] = 1
  sq1 = 1
  sq2 = 1
  q  = Array(Rational{compute}, n)
  m  = Array(Rational{compute}, n)
  b = 1
  while sq2 <= z
    # look for first non zero term in the series taylor(l)
    while l[b] == 0
      if (b += 1) > sl
        c = 1
        for k = 1:sq2
          c = lcm(c, den(q2[k]))
        end
        q = Array(output, sq2)
        q[:] = c*q2[1:sq2]
        return q
      end
    end
    # reciprocal of the series: x^(b-1) / taylor(l)
    m[1] = 1 // l[b]
    m[2:sl] = 0
    sm = 1
    for k = (b+1):sl
      c = compute(0)
      for j = b:(k-1)
        c -= l[j+1]*m[k-j]
      end
      m[ sm += 1 ] = c // l[b]
    end
    l, m = m, l
    sl = sm
    # compute (only) denominator of Pade approximant
    q[1:sq2] = l[1] * q2[1:sq2]
    q[sq2+1:b+sq1-1] = 0
    q[b:b+sq1-1] += q1[1:sq1]
    sq = max(sq2, b+sq1-1)
    while q[sq] == 0
      sq -= 1
    end
    q1, q2, q = q2, q, q1
    sq1, sq2 = sq2, sq
    # substract term of degree 0 before computing next reciprocal
    l[1] = 0
    b = 2
  end
  return Array(output, 0)

end

# println(recvec([1,1,2,3,5,8,13,21], compute=Int64, output=Int32))
# println(recvec([1,2,4,8,16,32,64]))
# println(recvec([64,32,16,8,4,2,1]))
# println("Failure: ",recvec([64,32,16,8,4,2,1,17]))
# println("Failure: ",recvec([64]))
