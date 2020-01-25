import math
import pygame
import os
import random


#sets the screen to the top left corner of the monitor
os.environ["SDL_VIDEO_WINDOW_POS"] = "%d, %d" % (5, 5)
pygame.init()

class Sphere:
        def __init__(self, name, radius, reflectivity, mass, elasticity, color, position):
            self.radius = radius
            self.mass = mass
            self.elasticity = elasticity
            self.color = color
            self.position = position
            self.velocity = (0, 0, 0)
            self.angle = math.pi/2
            self.type = "sphere"
            self.name = name
            self.numType = 2
            self.reflectivity = reflectivity

        def getName(self):
            return self.name

        def getRadius(self):
            return self.radius

        def getElasticity(self):
            return self.elasticity

        def setPosition(self, position):
            self.position = position

        def setColor(self, Color):
            self.color = Color

        def getColor(self):
            return self.color

        def move(self):
            self.position = (self.position[0] + self.velocity[0], self.position[1] + self.velocity[1], self.position[2] + self.velocity[2])


        def getNumType(self):
            return self.numType

        def getType(self):
            return self.type

        def getPosition(self):
            return self.position

        def setVelocity(self, v):
            self.velocity = v

        def getRelectivity(self):
            return self.reflectivity





# Window size, window area, camera location X, camera location y, camera location z, the length of the direction vectors, length of phi-hat vectors, the length of the theta-hat vecotrs, resolution, field of view on the x-y plane, field of view on the z-axis, degrees per "pixel" for the phi angle, degrees per "pixel" for the theta angle, speed controller for the camera, look sensitivity multiplier, fram rate target, pixel type(0 = circles, 1 = squares), adaptive FPS targeting on/off(true/false)
wS = 1100
wA = wS * wS
cameraLocX = 0
cameraLocY = 0
cameraLocZ = 0
r = 1
phiH = 1
thetaH = 1
phiDirec = 90
thetaDirec = 90
res = 1110
fovPHI = 90
fovTHETA = 90
dppPHI = fovPHI / math.sqrt(res)
dppTHETA = fovTHETA / math.sqrt(res)
speedMult = .1
turnMult = .1
frTarget = 10
pixlType = 1
adaptiveFPS = False

#The other one was getting long so here is the second grouping of starting variables
#music toggler
music = False

#light source coordinates
lx = 1
ly = 0
lz = 0

sphereLight = Sphere("LightSource", 0.1, 0, 100, 1, (255, 255, 255), (lx, ly, lz))



#point locations and colors
pX = [lx, 0, 1, -1]
pY = [ly, -1, 0, 0]
pZ = [lz, 0, 0, 0]
pColor = [(255, 255, 255), (0, 255, 0), (255, 0, 0), (0, 0, 255)]
#cube center points, side lengths, and colors
cubeX = [10]
cubeY = [10]
cubeZ = [10]
cubeS = [5]
cubeC = [(0, 0, 150)]

#sphere creations, (name, radii, reflectivity, mass, elasticity, color, position)
sphere1 = Sphere("worldObject", 2, 1, 100, 1, (127, 127, 127), (1, 2, 0))
sphere2 = Sphere("worldObject", 1, 1, 10, 1, (0, 0, 150), (3, 0, 0))
sphere3 = Sphere("worldObject", 0.5, 1, 5, 1, (255, 255, 255), (-1, 0, 0))


sphere1.setVelocity((0.00, 0.00, 0))


# all of the object lists
ob = [sphere1, sphere2, sphere3, sphereLight]

BLACK = (0, 0, 0)
WHITE = (255, 255, 255)


def sphereRunner(obj, camX, camY, camZ, phi, theta, oDist, oColor):
    position = ob[u].getPosition()
    if inView(cameraLocX, cameraLocY, position[0], position[1], curPHI):
        # boolean to check if the ray intersects the sphere, distance from camera to intersect point
        cS, distS, Intersect = checkSphere(cameraLocX, cameraLocY, cameraLocZ, ob[u], r, curPHI, curTHETA)
        # user can toggle the light source so that it can be like a "flashlight"
        if keys[pygame.K_g]:
            userL = True
        else:
            userL = False
        if cS and distS < oDist:
            sS = sphereShader(cameraLocX, cameraLocY, cameraLocZ, ob[u], r, curPHI, curTHETA, userL)
            color = sphereReflection(cameraLocX, cameraLocY, cameraLocZ, ob[u], Intersect, curPHI, curTHETA)
            if keys[pygame.K_t]:
                sS = 1
            oColor = (color[0] * sS, color[1] * sS, color[2] * sS)
            oDist = distS
        elif cS and oDist < 0:
            sS = sphereShader(cameraLocX, cameraLocY, cameraLocZ, ob[u], r, curPHI, curTHETA, userL)
            color = sphereReflection(cameraLocX, cameraLocY, cameraLocZ, ob[u], Intersect, curPHI, curTHETA)
            if keys[pygame.K_t]:
                sS = 1
            oColor = (color[0] * sS, color[1] * sS, color[2] * sS)

            oDist = distS
    return oColor, oDist

def cubeRunner(obj):
    pass

def pointRunner(obj):
    pass

def checkPoint(x0, y0, z0, x, y, z, r, phi, theta):
    phi = math.radians(phi)
    theta = math.radians(theta)

    a = (math.cos(phi) * math.sin(theta) * r) #+ (math.cos(phi) * math.cos(theta) * thetaH) - (math.sin(phi) * phiH)
    b = (math.sin(phi) * math.sin(theta) * r) #+ (math.sin(phi) * math.cos(theta) * thetaH) + (math.cos(phi) * phiH)
    c = (math.cos(theta) * r) #- (math.sin(theta) * thetaH)

    # x-combined with x0, y-combined with y0, z-combined with z0
    xC = x0 - x
    yC = y0 - y
    zC = z0 - z
#test comment
    #cross product the two matricies
    i = (yC * c) - (zC * b)
    j = (xC * c) - (zC * a)
    k = (xC * b) - (yC * a)

    d = math.sqrt(((i * i) + (j * j) + (k * k)) / ((a * a) + (b * b) + (c * c)))

    #distance from point to camera
    dd = math.sqrt(((x0 - x) * (x0 - x)) + ((y0 - y) * (y0 - y)) + ((z0 - z) * (z0 - z)))
    if d <= 0.05:
        return True, dd
    else:
        return False, dd

def checkCube():

    return False

def checkSphere(x0, y0, z0, sphere, r, phi, theta):
    phi = math.radians(phi)
    theta = math.radians(theta)

    pos = sphere.getPosition()
    h = pos[0]
    k = pos[1]
    l = pos[2]
    sR = sphere.getRadius()

    a = (math.cos(phi) * math.sin(theta) * r)  # + (math.cos(phi) * math.cos(theta) * thetaH) - (math.sin(phi) * phiH)
    b = (math.sin(phi) * math.sin(theta) * r)  # + (math.sin(phi) * math.cos(theta) * thetaH) + (math.cos(phi) * phiH)
    c = (math.cos(theta) * r)  # - (math.sin(theta) * thetaH)

    # a b and c for the quadratic involving the line and the sphere
    bb = ((2 * x0 * a) - (2 * h * a) + (2 * y0 * b) - (2 * k * b) + (2 * z0 * c) - (2 * l * c))
    aa = (a * a) + (b * b) + (c * c)
    cc = ((x0 * x0) + (y0 * y0) + (z0 * z0) + (h * h) + (k * k) + (l * l) - (2 * h * x0) - (2 * k * y0) - (2 * l * z0) - (sR * sR))

    # now do the formula of the discriminant (B^2 - 4AC)
    discrim = (bb * bb) - (4 * aa * cc)


    if discrim >= 0:



        # quadratic equation
        # its supposed to be -b + or - but I still need to implement that
        # actually I found that only the - seems to result in a sphere that is not inverted...this is strange as the - should only work for certain quadrants and angles
        t = ((-1 * bb) - math.sqrt((bb * bb) - (4 * aa * cc))) / (2 * aa)

        # now plug t in to the parametric forms of the line to get the intersect point(s)
        x = (t * a) + x0
        y = (t * b) + y0
        z = (t * c) + z0

        # find distance from camera or light source to intersect point
        d = math.sqrt(((x0 - x) * (x0 - x)) + ((y0 - y) * (y0 - y)) + ((z0 - z) * (z0 - z)))
        return True, d, (x, y, z)
    else:
        #placeholder distance (will never even be checked)
        d = 1000000
        return False, d, (0, 0, 0)

def sphereShader(x0, y0, z0, sphere, r, phi, theta, userLight):
    phi = math.radians(phi)
    theta = math.radians(theta)

    a = (math.cos(phi) * math.sin(theta) * r)  # + (math.cos(phi) * math.cos(theta) * thetaH) - (math.sin(phi) * phiH)
    b = (math.sin(phi) * math.sin(theta) * r)  # + (math.sin(phi) * math.cos(theta) * thetaH) + (math.cos(phi) * phiH)
    c = (math.cos(theta) * r)  # - (math.sin(theta) * thetaH)

    #we must now find t for the line equations to get the point(s) that they intersect at
    #a b and c for the quadratic involving the line and the sphere
    pos = sphere.getPosition()
    h = pos[0]
    k = pos[1]
    l = pos[2]
    sR = sphere.getRadius()

    bb = ((2 * x0 * a) - (2 * h * a) + (2 * y0 * b) - (2 * k * b) + (2 * z0 * c) - (2 * l * c))
    aa = (a * a) + (b * b) + (c * c)
    cc = ((x0 * x0) + (y0 * y0) + (z0 * z0) + (h * h) + (k * k) + (l * l) - (2 * h * x0) - (2 * k * y0) - (2 * l * z0) - (sR * sR))

    #quadratic equation
    #its supposed to be -b + or - but I still need to implement that
    #actually I found that only the - seems to result in a sphere that is not inverted
    t = ((-1 * bb) - math.sqrt((bb * bb) - (4 * aa * cc))) / (2 * aa)

    #now plug t in to the parametric forms of the line to get the intersect point(s)
    x = (t * a) + x0
    y = (t * b) + y0
    z = (t * c) + z0

    #the light source location
    llx = x0
    lly = y0
    llz = z0

    #find distance from camera or(/and) light source to intersect point
    d = math.sqrt(((llx - x) * (llx - x)) + ((lly - y) * (lly - y)) + ((llz - z) * (llz - z)))
    dd = math.sqrt(((lx - x) * (lx - x)) + ((ly - y) * (ly - y)) + ((lz - z) * (lz - z)))
    #checks to see if the light source is toggled to the user or not
    if userLight == False:
        d += dd
    #now return a fraction dependant on the size of the distance
    maxDist = 1

    distanceThreshold = 1
    if sphere.getName() == "LightSource":
        return (maxDist / (distanceThreshold))
    elif ((d*d) >= distanceThreshold):
        return (maxDist / (d * d))
    else:
        return (maxDist/ distanceThreshold)


def sphereReflection(x0, y0, z0, sphere, intersect, phi, theta):
    ref = sphere.getRelectivity()
    if ref == 0:
        return sphere.getColor()

    #see notebook for further explanations on the math
    phi = math.radians(phi)
    theta = math.radians(theta)

    #direction Vector of the original Ray
    a = (math.cos(phi) * math.sin(theta) * r)
    b = (math.sin(phi) * math.sin(theta) * r)
    c = (math.cos(theta) * r)

    #the intersect point
    xI = intersect[0]
    yI = intersect[1]
    zI = intersect[2]



    #center of the sphere
    pos = sphere.getPosition()
    h = pos[0]
    k = pos[1]
    l = pos[2]


    #the direction vector of the reflection line
    d = h - xI
    e = k - yI
    f = l - zI


    #Dot product of original vector and reflection vector
    dP = (d * a) + (e * b) + (f * c)

    #direction vector of the reflected line
    # http://www.3dkingdoms.com/weekly/weekly.php?a=2
    aR = -((2 * dP * d) - a) #xI - xR
    bR = -((2 * dP * e) - b) #yI - yR
    cR = -((2 * dP * f) - c) #zI - zR


    phiR = math.atan(bR/aR)

    col = sphere.getColor()

    for ur in range(0, len(ob)):
        #not ob[ur] == 0 is for points

        if not ob[ur] == sphere and not ob[ur] == 0:
            pos2 = ob[ur].getPosition()
            h2 = pos2[0]
            k2 = pos2[1]
            l2 = pos2[2]
            sR = ob[ur].getRadius()
            # a b and c for the quadratic involving the reflected line and other spheres
            bb = ((2 * xI * aR) - (2 * h2 * aR) + (2 * yI * bR) - (2 * k2 * bR) + (2 * zI * cR) - (2 * l2 * cR))
            aa = (aR * aR) + (bR * bR) + (cR * cR)
            cc = ((xI * xI) + (yI * yI) + (zI * zI) + (h2 * h2) + (k2 * k2) + (l2 * l2) - (2 * h2 * xI) - (
                    2 * k2 * yI) - (2 * l2 * zI) - (sR * sR))

            # now do the formula of the discriminant (B^2 - 4AC)
            discrim = (bb * bb) - (4 * aa * cc)

            colR = ob[ur].getColor()
            if discrim >= 0:
                if ob[ur].getName() == "LightSource":
                    return ob[ur].getColor()
                return (
                    ((colR[0] * ref) + col[0]) / 2, ((colR[1] * ref) + col[1]) / 2, ((colR[2] * ref) + col[2]) / 2)

    #returns color mixed with black void depending on its reflectivity
    return col #(col[0] * (1.1-ref), col[1] * (1.1-ref), col[2] * (1.1-ref))


#checks to make sure that the object is in front of the camera
def inView(x0, y0, x, y, phi):

    phiT = phi + 90
    phiT = math.radians(phiT)
    phi = math.radians(phi)

    s = math.tan(phiT)

    k = math.sin(phi)/math.fabs(math.sin(phi))

    return (k * y) >= k * ((s * (x - x0)) + y0)


window = pygame.display.set_mode((wS, wS), 0, 0)

# n is the diameter of each circle
n = math.sqrt(wA / res)
# pPerR is the pixels per row/column
pPR = wS / n

# setup mouse controls
pygame.mouse.set_pos(wS / 2, wS / 2)
pygame.event.set_grab(True)
pygame.mouse.set_visible(False)




font = pygame.font.Font(None, 30)
clock = pygame.time.Clock()

count = 0

if music == True:
    pygame.mixer.music.load("TheOuterWorlds.mp3")
    pygame.mixer.music.play()

while True:
    sphere1.move()
    if adaptiveFPS == True:
        if clock.get_fps() < frTarget:
            res -= 1
        else:
            res += 1

    dppPHI = fovPHI / math.sqrt(res)
    dppTHETA = fovTHETA / math.sqrt(res)
    # n is the diameter of each circle
    n = math.sqrt(wA / res)
    # pPerR is the pixels per row/column
    pPR = wS / n

    keys = pygame.key.get_pressed()
    window.fill(BLACK)

    # display fps
    fps = font.render(str(int(clock.get_fps())), True, WHITE)
    window.blit(fps, ((wS - 25), 5))
    clock.tick(60)

    # get the mouse movement vector values
    mous = pygame.mouse.get_rel()
    for event in pygame.event.get():
        pass
    if keys[pygame.K_ESCAPE]:
        break



    if keys[pygame.K_SPACE]:
        BLUE = (150, 150, 255)
        colorRand = (random.randint(0, 255), random.randint(0, 255), random.randint(0, 255))
        ob.append(Sphere("projectile", 1, 0.5, 5, 1, colorRand, (cameraLocX, cameraLocY, cameraLocZ)))




    if event.type == pygame.MOUSEMOTION:
        phiDirec -= (mous[0] * turnMult)
        thetaDirec += (mous[1] * turnMult)

    if keys[pygame.K_w]:
        cameraLocX = (math.cos(math.radians(phiDirec)) * math.sin(
            math.radians(thetaDirec)) * r * speedMult) + cameraLocX
        cameraLocY = (math.sin(math.radians(phiDirec)) * math.sin(
            math.radians(thetaDirec)) * r * speedMult) + cameraLocY
        cameraLocZ = (math.cos(math.radians(thetaDirec)) * r * speedMult) + cameraLocZ
    elif keys[pygame.K_s]:
        cameraLocX = (math.cos(math.radians(phiDirec)) * math.sin(
            math.radians(thetaDirec)) * -r * speedMult) + cameraLocX
        cameraLocY = (math.sin(math.radians(phiDirec)) * math.sin(
            math.radians(thetaDirec)) * -r * speedMult) + cameraLocY
        cameraLocZ = (math.cos(math.radians(thetaDirec)) * -r * speedMult) + cameraLocZ
    #for strafing add 90 degrees to the phi direction to go left and right, also we dont need the z movement because strafing only deals with the x-y plane
    if keys[pygame.K_a]:
        cameraLocX = (math.cos(math.radians(phiDirec + 90)) * math.sin(math.radians(thetaDirec)) * r * speedMult) + cameraLocX
        cameraLocY = (math.sin(math.radians(phiDirec + 90)) * math.sin(math.radians(thetaDirec)) * r * speedMult) + cameraLocY

    elif keys[pygame.K_d]:
        cameraLocX = (math.cos(math.radians(phiDirec + 90)) * math.sin(math.radians(thetaDirec)) * -r * speedMult) + cameraLocX
        cameraLocY = (math.sin(math.radians(phiDirec + 90)) * math.sin(math.radians(thetaDirec)) * -r * speedMult) + cameraLocY

    # I want the middle to be pointing straight so I added half the field of view to the starting phi angle
    curPHI = phiDirec + (fovPHI / 2)
    # I want the middle to be pointing straight so I subtracted half the theta field of view to make it start facing up
    curTHETA = thetaDirec - (fovTHETA / 2)

    #for circles
    if pixlType == 0:
        curX = n / 2
        curY = n / 2
    elif pixlType == 1:
        curX = 0
        curY = 0
    #for squares


    for f in range(int(pPR)):
        for l in range(int(pPR)):

            #distance from camera to closest object that has intersect point with camera
            objDist = -1
            #Color of closest object to camera that intersects with ray
            objC = (0, 0, 0)
            for u in range(0, len(ob)):
                numType = ob[u].getNumType()
                types = [pointRunner, cubeRunner, sphereRunner]
                objC, objDist = types[numType](ob[u], cameraLocX, cameraLocY, cameraLocZ, curPHI, curTHETA, objDist, objC)



            #draw the pixel based on which object is closest
            if objDist >= 0:

                if pixlType == 0:
                    pygame.draw.circle(window, objC, (int(curX), int(curY)), int(n / 2))
                elif pixlType == 1:
                    pixlRect = pygame.Rect(curX, curY, n, n)
                    pygame.draw.rect(window, objC, pixlRect)


            curPHI -= dppPHI
            curX += n
        #for circles
        if pixlType == 0:
            curX = n / 2
        elif pixlType == 1:
            curX = 0
        curY += n

        curPHI = phiDirec + (fovPHI / 2)
        curTHETA += dppTHETA

    w = font.render("W", True, (150, 150, 150))
    a = font.render("A", True, (150, 150, 150))
    s = font.render("S", True, (150, 150, 150))
    d = font.render("D", True, (150, 150, 150))
    if keys[pygame.K_w]:
        w = font.render("W", True, WHITE)

    if keys[pygame.K_a]:
        a = font.render("A", True, WHITE)

    if keys[pygame.K_s]:
        s = font.render("S", True, WHITE)

    if keys[pygame.K_d]:
        d = font.render("D", True, WHITE)

    window.blit(w, (30, (wS - 50)))
    window.blit(a, (5, (wS - 25)))
    window.blit(s, (30, (wS - 25)))
    window.blit(d, (55, (wS - 25)))

    pygame.display.update()

pygame.quit()
quit()
