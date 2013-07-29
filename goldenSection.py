import math
phi = 0.5 * (1 + math.sqrt(5));
resphi = 2 - phi;

def goldenSectionSearch(f, p, epsilon):
    print p
    a,b,c = zip(*p)[0]
    x = ( (b + resphi * (c - b))  if (c - b > b - a) else
          (b - resphi * (b - a)) )
    if abs(c-a) < epsilon:
        return 0.5 * (c + a)
    fb = p[1][1]
    if fb==None: fb = f(b)
    fx = f(x)
    assert fx != fb
    if fx < fb:
        return ( goldenSectionSearch(f, (p[1], (x,fx), p[2]), epsilon) if (c - b > b - a) else
                 goldenSectionSearch(f, (p[0], (x,fx), p[1]), epsilon) )
    return ( goldenSectionSearch(f, (p[0], p[1], (x,fx)), epsilon) if (c - b > b - a) else 
             goldenSectionSearch(f, ((x,fx), p[1], p[2]), epsilon) )
