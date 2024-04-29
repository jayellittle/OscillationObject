// oscillation7

// physical parameters

float x0 = 0.05;
float y0 = 0.2;

float m = 0.5;

float radius = 0.05;

float k = 10000.0;
float l0 = 0.2;
float d = 100.0;

float wallY = 3.0;
float wallK = 10000.0;
float wallD = 100.0;

float gx = 0.0;
float gy = 9.8;

float dt = 3e-6;

// viewing parameters

float viewingSize = 4.0;

int xOffset;
int yOffset;
float viewingScale;

// objects

int nParticles = 5;
int nSprings = nParticles - 1;

Particle[] particles;
Spring[] springs;
Wall wall;

float gridSizeX = 0.2;
float gridSizeY = 0.2;

int gridNumX = 8;
int gridNumY = 8;

float gridOffsetX;
float gridOffsetY;

// particle //////////////////////////////

class Particle {
  float x, y;
  float vx, vy;
  float fx, fy;
  
  float m;
  boolean isFixed;
  
  float radius;
  
  float xPrev, yPrev;
  float vxPrev, vyPrev;
  float fxPrev, fyPrev;
  
  Particle(float x0, float y0, float vx0, float vy0, float m0,
           boolean f, float r) {
    m = m0;
    x = x0;
    y = y0;
    vx = vx0;
    vy = vy0;
    isFixed = f;
    radius = r;
  }
  
  void init() {
    fx = 0.0;
    fy = 0.0;
    
    xPrev = x;
    yPrev = y;
    vxPrev = vx;
    vyPrev = vy;
    fxPrev = fx;
    fyPrev = fy;
  }
  
  void clearForce() {
    fx = 0.0;
    fy = 0.0;
  }
  
  void addForce(float x, float y) {
    fx += x;
    fy += y;
  }
  
  void move(float dt) {
    if (isFixed) {
      return;
    }
    
    float xNew = x + (3 * vx - vxPrev) * 0.5 * dt;
    float yNew = y + (3 * vy - vyPrev) * 0.5 * dt;
    float vxNew = vx + (3 * fx - fxPrev) * 0.5 * dt / m;
    float vyNew = vy + (3 * fy - fyPrev) * 0.5 * dt / m;
    
    xPrev = x;
    yPrev = y;
    vxPrev = vx;
    vyPrev = vy;
    fxPrev = fx;
    fyPrev = fy;
    
    x = xNew;
    y = yNew;
    vx = vxNew;
    vy = vyNew;
  }
  
  void draw() {
    if (radius <= 0.0) {
      return;
    }
    
    if (isGrabbed == true && this == grabbedParticle) {
      fill(240, 128, 128);
    }
    
    fill(224);
    stroke(128);
    strokeWeightScaled(1.0);
    
    pushMatrix();
    translate(x, y);
    ellipse(0, 0, radius, radius);
    popMatrix();
  }
};

// spring ////////////////////////////////////

class Spring {
  Particle[] particles;
  float k;
  float l0;
  float d;
  
  Spring(Particle p0, Particle p1, float k0, float l00, float d0)
  {
    particles = new Particle[2];
    particles[0] = p0;
    particles[1] = p1;
    k = k0;
    l0 = l00;
    d = d0;
  }
  
  void init() {
  }
  
  void calc() {
    float dx = particles[1].x - particles[0].x;
    float dy = particles[1].y - particles[0].y;
    float l = sqrt(dx * dx + dy * dy);
    float ex = dx / l;
    float ey = dy / l;
    float fx = -k * (l - l0) * ex;
    float fy = -k * (l - l0) * ey;
    
    float vx = particles[1].vx - particles[0].vx;
    float vy = particles[1].vy - particles[0].vy;
    fx += -d * vx;
    fy += -d * vy;
    
    particles[0].addForce(-fx, -fy);
    particles[1].addForce(fx, fy);
  }
  
  void draw() {
    noFill();
    stroke(128);
    strokeWeightScaled(1.0);
    
    line(particles[0].x, particles[0].y, particles[1].x, particles[1].y);
  }
};

// wall ///////////////////////////////////////////

class Wall {
  float y;
  float k;
  float d;
  
  Wall(float y0, float k0, float d0) {
    y = y0;
    k = k0;
    d = d0;
  }
  
  float collisionForceX(float dy, float dvx, float dvy) {
    return -d * dvx;
  }
  
  float collisionForceY(float dy, float dvx, float dvy) {
    return -k * dy - d * dvy;
  }
  
  void init() {
  }
  
  void calc() {
  }
  
  void draw() {
    fill(224);
    noStroke();
    rect(-viewingSize / 2, y, viewingSize, viewingSize);
    
    stroke(128);
    strokeWeightScaled(1.0);
    line(-viewingSize, y, viewingSize, y);
  }
};

// simulation /////////////////////////////////////

void simulationInit() {
  gridOffsetX = -gridSizeX * (gridNumX - 1) / 2.0;
  gridOffsetY = -gridSizeY * (gridNumY - 1) / 2.0;
  
  nParticles = gridNumX * gridNumY;
  nSprings = gridNumX * (gridNumY - 1) + (gridNumX - 1) * gridNumY
      + 2 * (gridNumX - 1) * (gridNumY - 1);

  particles = new Particle[nParticles];
  for (int j = 0; j < gridNumY; j++) {
    for (int i = 0; i < gridNumX; i++) {
      int n = i + j * gridNumX;
      particles[n] = new Particle(gridOffsetX + i * gridSizeX,
                                  gridOffsetY + j * gridSizeY,
                                  0.0, 0.0, m, false, radius);
    }
  }
  
  springs = new Spring[nSprings];
  int n = 0;
  for (int j = 0; j < gridNumY; j++) {
    for (int i = 0; i < gridNumX - 1; i++) {
      springs[n] = new Spring(particles[i + j * gridNumX],
                              particles[(i + 1) + j * gridNumX],
                              k, gridSizeX, d);
      n++;
    }
  }
  
  for (int j = 0; j < gridNumY - 1; j++) {
    for (int i = 0; i < gridNumX; i++) {
      springs[n] = new Spring(particles[i + j * gridNumX],
                              particles[i + (j + 1) * gridNumX],
                              k, gridSizeY, d);
      n++;
    }
  }
  
  float l = sqrt(gridSizeX * gridSizeX + gridSizeY * gridSizeY);
  for (int j = 0; j < gridNumY - 1; j++) {
    for (int i = 0; i < gridNumX - 1; i++) {
      springs[n] = new Spring(particles[i + j * gridNumX],
                              particles[(i + 1) + (j + 1) * gridNumX],
                              k, l, d);
      n++;
    }
  }
  for (int j = 0; j < gridNumY - 1; j++) {
    for (int i = 0; i < gridNumX - 1; i++) {
      springs[n] = new Spring(particles[(i + 1) + j * gridNumX],
                              particles[i + (j + 1) * gridNumX],
                              k, l, d);
      n++;
    }
  }
  
  for (Particle p: particles) {
    p.init();
  }
  
  for (Spring s: springs) {
    s.init();
  }
  
  wall = new Wall(wallY, wallK, wallD);
  wall.init();
}

void simulationCalc() {
  while (true) {
    if (isGrabbed) {
      float mx = (mouseX - xOffset) / viewingScale;
      float my = (mouseY - yOffset) / viewingScale;
      grabbedParticle.x = mx;
      grabbedParticle.y = my;
      
      float vx = (mouseX - pmouseX) / viewingScale;
      float vy = (mouseY - pmouseY) / viewingScale;
      grabbedParticle.vx = vx;
      grabbedParticle.vy = vy;
    }
    
    for (Particle p: particles) {
      p.clearForce();
      p.addForce(p.m * gx, p.m * gy);
    }
    
    for (Particle p: particles) {
      float dy = (p.y + p.radius) - wall.y;
      if (dy > 0.0) {
        p.addForce(wall.collisionForceX(dy, p.vx, p.vy),
                   wall.collisionForceY(dy, p.vx, p.vy));
      }
    }
    
    for (Spring s: springs) {
      s.calc();
    }
    
    for (Particle p: particles) {
      p.move(dt);
    }
  }
}

void simulationDraw() {
  ellipseMode(RADIUS);
  background(255);
  
  translate(xOffset, yOffset);
  scale(viewingScale);
  
  wall.draw();
  
  for (Spring s: springs) {
    s.draw();
  }
  
  for (Particle p: particles) {
    p.draw();
  }
}

// utility function ///////////////////////////////

void strokeWeightScaled(float s) {
  strokeWeight(s / viewingScale);
}

// setup //////////////////////////////////////////

void setup() {
  size(512, 512);
  frameRate(60);
  smooth(4);
  
  xOffset = width / 2;
  yOffset = 0;
  viewingScale = width / viewingSize;
  
  simulationInit();
  
  thread("simulationCalc");
}

// draw ///////////////////////////////////////////

void draw() {
  simulationDraw();
}

// callbacks //////////////////////////////////////

void keyPressed() {
  if (key == ESC || key == 'q') {
    exit();
  }
}

// interaction ////////////////////////////////////

boolean isGrabbed = false;
Particle grabbedParticle;
float grabStartX;
float grabStartY;

void mousePressed() {
  float mx = (mouseX - xOffset) / viewingScale;
  float my = (mouseY - yOffset) / viewingScale;
  for (Particle p: particles) {
    float dx = p.x - mx;
    float dy = p.y - my;
    if (dx * dx + dy * dy < p.radius * p.radius) {
      isGrabbed = true;
      grabbedParticle = p;
      grabbedParticle.isFixed = true;
      grabStartX = mouseX;
      grabStartY = mouseY;
      break;
    }
  }
}

void mouseReleased() {
  if (isGrabbed) {
    isGrabbed = false;
    grabbedParticle.isFixed = false;
  }
}
