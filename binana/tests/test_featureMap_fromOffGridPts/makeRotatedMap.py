import numpy
import peel

def normalize(v, tolerance=0.00001):
    mag2 = sum(n * n for n in v)
    if abs(mag2 - 1.0) > tolerance:
        mag = numpy.sqrt(mag2)
        v = tuple(n / mag for n in v)
    return v

def q_mult(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
    z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
    return w, x, y, z

def q_conjugate(q):
    
    q = normalize(q)
    w, x, y, z = q
    return (w, -x, -y, -z)


def qv_mult(q1, v1):
    q2 = (0.0,) + v1
    return q_mult(q_mult(q1, q2), q_conjugate(q1))[1:]






origFrame = numpy.load('POVME_frame_1.npy')
origFrame = origFrame - [int(numpy.mean(origFrame[:,0])),
                         int(numpy.mean(origFrame[:,1])),
                         int(numpy.mean(origFrame[:,2]))]
#print origFrame

rotQuad =  (0.81634468808940908, -0.45485698954637466, -0.33032253549642465, 0.13256504755110923)

rotFrame = []
for coord in origFrame:
    coord = tuple(coord)
    coord = list(qv_mult(rotQuad,coord))
    rotFrame.append(coord)
rotFrame = numpy.array(rotFrame)

#print rotFrame
origFrameFM = peel.featureMap.fromPovmeList(origFrame, 1, justCoords=True)
#rotFrameFM = peel.featureMap.fromPovmeList(rotFrame, 1, justCoords=True)
for i in range(5):
    rotFrame_FOGP_FM = peel.featureMap.fromOffGridPts(rotFrame, 1, justCoords=True)

origFrameFM.write_pdb('origFrame.pdb')
with open('rotFrame.pdb','wb') as of:
    of.write(peel.numpy_to_pdb(rotFrame,'X'))
#rotFrameFM.write_pdb('rotFrame.pdb')
rotFrame_FOGP_FM.write_pdb('rotFrame_FOGP.pdb')

