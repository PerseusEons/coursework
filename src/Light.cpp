#include <CanvasPoint.h>
#include <CanvasTriangle.h>
#include <Colour.h>
#include <DrawingWindow.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>
#include <TextureMap.h>
#include <Utils.h>

#include <fstream>
#include <functional>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;
using namespace glm;
#define max(x,y) ((x > y) ? (x) : y)
#define min(x,y) ((x < y) ? (x) : y)
#define pb push_back
#define abs(x) (((x) > 0) ? (x) : -(x))
#define double float
#define WIDTH 600
#define tri triangle
#define vert vertices
#define col colour
#define CP CanvasPoint
#define HEIGHT 600
#define ttP texturePoint
const double PI = acosl(-1);
const Colour base_col = Colour(255, 255, 255);


vec3 o(0.0, 0.0, 0.0);
vec3 cam(0.0, 0.0, 4.0);
vec3 light(0.0, 1.0, 2.0);
mat3 cam_orientation(vec3(1.0, 0.0, 0.0), 
                    vec3(0.0, 1.0, 0.0),
                    vec3(0.0, 0.0, 1.0));
#define wsPIC window.setPixelColour
double foc = 500.0;
bool orbiting = false, show_light = false;
bool proximity = true, angle_of = true, shadows = true, specular = true;

void draw(DrawingWindow &window) {window.clearPixels(); }

void update(DrawingWindow &window) {
  // Function for performing animation (shifting artifacts or moving the camera)
}

vector<CP> inter_points_in_segment(CP start, CP to,
                                            int step) {
  vector<CP> points = {start};

  CP tmp = start;
  for (int i = 0; i < step - 1; i++) {
    double sx = start.x;
    double sy = start.y;
    tmp.x = tmp.x + (to.x - sx) / (step - 1);
    tmp.y = tmp.y + (to.y - sy) / (step - 1);
    tmp.depth = tmp.depth + (to.depth - start.depth) / (step - 1);
    points.pb(tmp);
  }
  return points;
}

vector<TexturePoint> inter_points_in_segment(TexturePoint start,
                                             TexturePoint to, int step) {
  vector<TexturePoint> points = {start};

  TexturePoint tmp = start;
  for (int i = 0; i < step - 1; i++) {
    double sx = start.x;
    double sy = start.y;
    tmp.x = tmp.x + (to.x - sx) / (step - 1);
    tmp.y = tmp.y + (to.y - sy) / (step - 1);
    points.pb(tmp);
  }
  return points;
}
bool ok(double x, double y, double w, double h) {
  if (x >= 0 && x < w && y >= 0 && y < h) {
    return 1;
  }
  return 0;
}
void d_l(CP start, CP to, Colour col,
               DrawingWindow &window) {
  double ww = window.width;
  double wh = window.height;
  double sx = start.x;
  double sy = start.y;
  double x_dis = to.x - sx;
  double y_dis = to.y - sy;
  double step = max(abs(x_dis), abs(y_dis));

  uint32_t c = (255 <<24) + (int(col.red) * 256 * 256) +
               (int(col.green) * 256) + int(col.blue);

  for (double i = 0.0; i < step; i++) {
    int x = roundl(sx + (x_dis / step * i));
    int y = roundl(sy + (y_dis / step * i));
    if (ok(x, y, ww, wh)) {
      wsPIC(x, y, c);
    }
  }
}

void draw_tri(CanvasTriangle tri, Colour col,
                   DrawingWindow &window) {

  d_l(tri[0], tri[1], col, window);
  d_l(tri[1], tri[2], col, window);
  d_l(tri[2], tri[0], col, window);
}

CP find_bot1(CP top, CP bot1, CP bot2) {
  double b1y = bot1.y;
  double b2y = bot2.y;
  double b2x = bot2.x;
  double tpy = top.y;
  double tpx = top.x;
  double k = ((b1y - tpy) / (b2y - tpy)) * (b2x - tpx);
  double x_bot1 = tpx + k;
  return CP(roundl(x_bot1), b1y);
}

void f_h_tri(CanvasTriangle tri, Colour col, DrawingWindow &window,
             vector<vector<double>> &depths) {
  double ww = window.width;
  double wh = window.height;
  CP top = tri.vert[0];
  CP bot1 = tri.vert[1];
  CP bot2 = tri.vert[2];

  double b1y = bot1.y;
  double b2y = bot2.y;
  double tpy = top.y;
  vector<CP> left =
      inter_points_in_segment(top, bot1, abs(b1y - tpy) + 2);
  vector<CP> right =


      inter_points_in_segment(top, bot2, abs(b2y - tpy) + 2);

  for (int i = 0; i < left.size(); i++) {
    int step = abs(left[i].x - right[i].x);

    vector<CP> points =
        inter_points_in_segment(left[i], right[i], step + 2);
    for (auto j : points) {
      int x = roundl(j.x);
      int y = roundl(j.y);
      if (ok(x, y, ww, wh)) {
        if (-1 / j.depth > depths[x][y]) {
          depths[x][y] = -1 / j.depth;
          uint32_t c = (255 <<24) + (int(col.red) * 256 * 256) +
                       (int(col.green) * 256) + int(col.blue);
          wsPIC(x, y, c);
        }
      }
    }
  }
}

void texture_half_tri(CanvasTriangle tri, TextureMap texture,
                           DrawingWindow &window,


                           vector<vector<double>> &depths) {
  CP top = tri.vert[0];
  CP bot1 = tri.vert[1];
  CP bot2 = tri.vert[2];
  double b1y = bot1.y;
  double b2y = bot2.y;
  double tpy = top.y;
  vector<CP> left =
      inter_points_in_segment(top, bot1, abs(b1y - tpy) + 2);
  vector<CP> right =
      inter_points_in_segment(top, bot2, abs(b2y - tpy) + 2);
  vector<TexturePoint> left_texture = inter_points_in_segment(
      top.ttP, bot1.ttP, abs(b1y - tpy) + 2);
  vector<TexturePoint> right_texture = inter_points_in_segment(


      top.ttP, bot2.ttP, abs(b1y - tpy) + 2);
  int ii=-1;
  for (auto i:left) {
    ii += 1;
    int step = abs(i.x - right[ii].x);

    vector<CP> points =
        inter_points_in_segment(i, right[ii], step + 2);
    vector<TexturePoint> points_texture =
        inter_points_in_segment(left_texture[ii], right_texture[ii], step + 2);
    int jj = -1;
    for (auto j : points) {
      jj += 1;
      int x = roundl(j.x);
      int y = roundl(j.y);
      if (x >= 0 && x < window.width && y >= 0 && y < window.height) {
        if (-1 / j.depth > depths[x][y]) {
          depths[x][y] = -1 / j.depth;


          wsPIC(x, y,
                texture.pixels[roundl(points_texture[jj].y) * texture.width +
                               roundl(points_texture[jj].x)]);
        }
      }
    }
  }
}

void f_tri(CanvasTriangle tri, Colour col, DrawingWindow &window,
           vector<vector<double>> &depths) {
  CP top = tri.vert[0];
  CP bot1 = tri.vert[1];
  CP bot2 = tri.vert[2];

  double b1y = bot1.y;
  double b2y = bot2.y;
  double tpy = top.y;
  if (b2y < b1y) swap(bot2, bot1);
  if (b1y < tpy) swap(bot1, top);
  if (b2y < b1y) swap(bot2, bot1);
  b1y = bot1.y;
  b2y = bot2.y;
  tpy = top.y;
  CP bot1_2 = find_bot1(top, bot1, bot2);
  bot1_2.depth = top.depth + (bot2.depth - top.depth) *
                                 ((b1y - tpy) / (b2y - tpy));



  CanvasTriangle tri_1 = CanvasTriangle(top, bot1, bot1_2);
  CanvasTriangle tri_2 = CanvasTriangle(bot2, bot1, bot1_2);

  f_h_tri(tri_1, col, window, depths);
  f_h_tri(tri_2, col, window, depths);

  // draw_tri(tri, base_col, window);
}

void texture_tri(TextureMap texture, CanvasTriangle tri,
                      DrawingWindow &window, vector<vector<double>> &depths) {
  CP top = tri.vert[0];
  CP bot1 = tri.vert[1];
  CP bot2 = tri.vert[2];

  double b1y = bot1.y;
  double b2y = bot2.y;
  double tpy = top.y;
  if (b2y < b1y) swap(bot2, bot1);
  if (b1y < tpy) swap(bot1, top);
  if (b2y < b1y) swap(bot2, bot1);
  b1y = bot1.y;
  b2y = bot2.y;
  tpy = top.y;

  CP bot1_2 = find_bot1(top, bot1, bot2);
  double ff = (b1y - tpy) ;
  double fff = (b2y - tpy) ;
  double btd = (bot2.depth - top.depth);
  bot1_2.depth = top.depth + (ff/ fff) * btd;

  bot1_2.ttP.x =
      top.ttP.x + (bot2.ttP.x - top.ttP.x) *
                               (b1y - tpy) / (b2y - tpy);
  bot1_2.ttP.y =
      top.ttP.y + (bot2.ttP.y - top.ttP.y) *
                               (by - tpy) / (b2y - tpy);

  texture_half_tri(CanvasTriangle(top, bot1, bot1_2), texture, window, depths);
  texture_half_tri(CanvasTriangle(bot2, bot1, bot1_2), texture, window, depths);

}

void draw_wireframe(vector<ModelTriangle> tris, DrawingWindow &window) {
  double ww = window.width;
  double wh = window.height;
  for (auto i : tris) {
    ModelTriangle tri = i;
    CanvasTriangle t;
    for (int j = 0; j < tri.vert.size(); j++) {
      vec3 ve = tri.vert[j];

      double xxx = ve.x - cam.x;
      double yyy = ve.y - cam.y;
      double zzz = ve.z - cam.z;
      vec3 cam_to_vertex(xxx, yyy, zzz);
      vec3 adj_v = cam_to_vertex * cam_orientation;
      double xx = (foc * (adj_v.x) / (adj_v.z));
      double yy = (foc * (adj_v.y) / (adj_v.z));
      int u = -xx + (ww / 2);
      int v = yy + (wh / 2);


      t.vert[j] = CP(u, v, adj_v.z);
      t.vert[j].ttP = tri.texturePoints[j];
    }

    draw_tri(t, base_col, window);
  }
  double xxx=light.x - cam.x,yyy=light.y - cam.y,zzz=light.z - cam.z;
  vec3 cam_to_vertex = vec3(xxx, yyy, zzz);
  vec3 adj_v = cam_to_vertex * cam_orientation;

  int xx = -(foc * (adj_v.x) / (adj_v.z)) + (ww / 2);
  int yy = (foc * (adj_v.y) / (adj_v.z)) + (wh / 2);

  // prints red PIxels to show light location
  wsPIC(xx, yy, (255 <<24) + (255 * 256 * 256) + (0 * 256) + 0);
  wsPIC(xx - 1, yy, (255 <<24) + (255 * 256 * 256) + (0 * 256) + 0);
  wsPIC(xx, yy - 1, (255 <<24) + (255 * 256 * 256) + (0 * 256) + 0);
  wsPIC(xx + 1, yy, (255 <<24) + (255 * 256 * 256) + (0 * 256) + 0);
  wsPIC(xx, yy + 1, (255 <<24) + (255 * 256 * 256) + (0 * 256) + 0);
}

void draw_rasterise(vector<ModelTriangle> tris, DrawingWindow &window) {
  double ww = window.width;
  double wh = window.height;
  vector<vector<double>> depths(
      ww, vector<double>(wh, -(numeric_limits<double>::infinity())));

  for (auto i : tris) {
    ModelTriangle tri = i;
    CanvasTriangle t;
    for (int j = 0; j < tri.vert.size(); j++) {
      vec3 ve = tri.vert[j];
      double xx = ve.x - cam.x;
      double yy = ve.y - cam.y;
      double zz = ve.z - cam.z;
      vec3 cam_to_vertex(xx, yy, zz);
      vec3 adj_v = cam_to_vertex * cam_orientation;

      t.vert[j] =
          CP(-(foc * (adj_v.x) / (adj_v.z)) + (ww / 2.0),


                      (foc * (adj_v.y) / (adj_v.z)) + (wh / 2.0), adj_v.z);
      t.vert[j].ttP = tri.texturePoints[j];
    }
    if (tri.col.name != "") {
      TextureMap texture(tri.col.name);
      for (int j = 0; j < t.vert.size(); j++) {
        t.vert[j].ttP.x =
            t.vert[j].ttP.x * texture.width;
        t.vert[j].ttP.y =


            t.vert[j].ttP.y * texture.height;
      }
      texture_tri(texture, t, window, depths);
    } else
      f_tri(t, tri.col, window, depths);
  }
}
bool in_range(double x,double l=0.0,double r=1.0) {
  if(x>=l && x<=r) return 1;
  return 0;
}
bool is_sha(RayTriangleIntersection intersect, vector<ModelTriangle> tris) {
  vec3 shadow_ray = light - intersect.intersectionPoint;
  int ii=-1;
  for (auto i : tris) {
    ii += 1;
    ModelTriangle tri = i;



    vec3 e0 = tri.vert[1] - (*tri.vert.begin());
    vec3 e1 = tri.vert[2] - (*tri.vert.begin());
    vec3 sp_vector = intersect.intersectionPoint - (*tri.vert.begin());
    mat3 de_matrix(-normalize(shadow_ray), e0, e1);
    vec3 possible_s = inverse(de_matrix) * sp_vector;
    double t = possible_s.x, u = possible_s.y, v = possible_s.z;

    if (in_range(u) && in_range(v) && (u + v) <= 1.0)
      if (t < glm::length(shadow_ray) && t > 0.05 &&
          ii != intersect.triangleIndex)
        return true;
  }
  return false;
}

double get_scale(RayTriangleIntersection rt_int, int scale) {
  vec3 normal = normalize(rt_int.intersectedTriangle.normal);
  vec3 light_ray = light - rt_int.intersectionPoint;
  vec3 view_ray = normalize(cam - rt_int.intersectionPoint);
  vec3 reflection_ray =
      normalize(normalize(light_ray) -
                (normal * 2.0f * dot(normalize(light_ray), normal)));

  double scale_p =
      (proximity) ? 40.0 / (4.0 * PI * (powl(length(light_ray), 2.0))) : 0.0;
  double s_a = dot(normal, normalize(light_ray));
  double s_s = pow(dot(reflection_ray, view_ray), scale);

  if (angle_of && s_a > 0 ) scale_p *= s_a;
  if (specular && s_s > 0 ) scale_p += s_s;
  return (scale_p < 1) ? scale_p : 1;
}



double gourad(RayTriangleIntersection rt_int, int scale) {
  ModelTriangle t = rt_int.intersectedTriangle;
  vec3 light_ray = light - rt_int.intersectionPoint;
  vec3 view_ray = normalize(cam - rt_int.intersectionPoint);

  vector<double> scales;
  vector<vec3> reflections;
  for (auto i:t.normals) {
    vec3 reflection_ray = normalize(
        normalize(light_ray) -
        (i * 2.0f * dot(normalize(light_ray), i)));
    reflections.pb(reflection_ray);

    double tmp_a = (angle_of) ? dot(i, normalize(light_ray)) : 1;
    scales.pb(tmp_a);
  }
  vec3 reflection_ray = (1 - rt_int.u - rt_int.v) * reflections[0] +
                        rt_int.u * reflections[1] + rt_int.v * reflections[2];

  double s_s = pow(dot(normalize(reflection_ray), view_ray), scale);


  double s_a = (1 - rt_int.u - rt_int.v) * scales[0] +
                   rt_int.u * scales[1] + rt_int.v * scales[2];

  double scale_p =
      (proximity) ? 40.0 * s_a / (4.0 * PI * (powl(length(light_ray), 2.0))) : 0.0;
  if (s_s > 0 && specular) scale_p += (s_s * 0.2);

  return (scale_p < 1) ? scale_p : 1;
}



double phong(RayTriangleIntersection rt_int, int scale) {
  ModelTriangle t = rt_int.intersectedTriangle;
  vec3 normal = (1 - rt_int.u - rt_int.v) * t.normals[0] +
                rt_int.u * t.normals[1] + rt_int.v * t.normals[2];
  vec3 light_ray = light - rt_int.intersectionPoint;
  vec3 view_ray = normalize(cam - rt_int.intersectionPoint);
  vec3 reflection_ray =
      normalize(normalize(light_ray) -
                (normal * 2.0f * dot(normalize(light_ray), normal)));

  double s_a = (angle_of) ? dot(normal, normalize(light_ray)) : 1;
  double scale_p =
      (proximity) ? 40.0 * s_a / (4.0 * PI * (powl(length(light_ray), 2.0))) : 0.0;


  double s_s = powl(dot(reflection_ray, view_ray), scale);
  if (s_s > 0 && specular) scale_p += (s_s * 0.2);
  return (scale_p < 1) ? scale_p : 1;
}

RayTriangleIntersection get_closest_intersection(vec3 direction,
                                                 vector<ModelTriangle> tris) {
  RayTriangleIntersection rti;
  rti.distanceFromCamera = numeric_limits<double>::infinity();
  vec3 ray = cam - direction;
  ray = normalize(cam_orientation * ray);
  int ii = -1;
  for (auto i :tris) {
    ii += 1;
    ModelTriangle tri = i;

    vec3 e0 = tri.vert[1] - tri.vert[0];
    vec3 e1 = tri.vert[2] - tri.vert[0];
    vec3 sp_vector = cam - tri.vert[0];
    mat3 de_matrix(-ray, e0, e1);
    vec3 possible_s = inverse(de_matrix) * sp_vector;
    double t = possible_s.x, u = possible_s.y, v = possible_s.z;



    if ((u + v) <= 1.0 && (u >= 0.0) && (u <= 1.0) && (v >= 0.0) && (v <= 1.0) )
      if (rti.distanceFromCamera > t && t > 0) {
        rti.distanceFromCamera = t;
        rti.intersectedTriangle = tri;
        // rti.intersectedTriangle.normal = cross(e1,e0);
        rti.triangleIndex = ii;
        rti.u = u;
        rti.v = v;

        vec3 intersect = tri.vert[0] + u * e0 + v * e1;
        rti.intersectionPoint = intersect;
      }
  }
  return rti;
}

function<double(RayTriangleIntersection rt_int, int scale)> brightness = phong;

void draw_raytrace(vector<ModelTriangle> tris, DrawingWindow &window) {
  double ww = window.width;
  double wh = window.height;
  for (int x = 0; x < ww; x++) {
    for (int y = 0; y < wh; y++) {
      RayTriangleIntersection rt_int = get_closest_intersection(
          vec3((int(ww) / 2) - x, y - (int(wh) / 2), foc), tris);



      double scale = brightness(rt_int, 64);
      scale = (scale > 0.15) ? scale : 0.15;
      if (!isinf(rt_int.distanceFromCamera)) {
        Colour col = rt_int.intersectedTriangle.col;
        uint32_t c = (255 <<24) + (int(col.red * scale) * 256 * 256) +
                     (int(col.green * scale) * 256 * 256) +
                     int(col.blue * scale);

        if (is_sha(rt_int, tris) && shadows) {
          double s_s = 0.1;
          uint32_t s = (255 <<24) + (int(col.red * s_s) * 256 * 256) +
                       (int(col.green * s_s) * 256 * 256) +
                       int(col.blue * s_s);
          wsPIC(x, y, s);
        } else
          wsPIC(x, y, c);
      }
    }
  }
  if (show_light) {
    double lx = light.x,ly = light.y,lz = light.z;
    double cx = cam.x,cy = cam.y,cz = cam.z;
    vec3 cam_to_vertex =
        vec3(lx - cx, ly - cy, lz - cz);
    vec3 adj_v = cam_to_vertex * cam_orientation;

    int xx = -(foc * (adj_v.x) / (adj_v.z)) + (ww / 2);
    int yy = (foc * (adj_v.y) / (adj_v.z)) + (wh / 2);

    wsPIC(xx, yy, (255 <<24) + (255 * 256 * 256) + (0 * 256 * 256) + 0);
    wsPIC(xx + 1, yy, (255 <<24) + (255 * 256 * 256) + (0 * 256 * 256) + 0);
    wsPIC(xx, yy + 1, (255 <<24) + (255 * 256 * 256) + (0 * 256 * 256) + 0);
    wsPIC(xx - 1, yy, (255 <<24) + (255 * 256 * 256) + (0 * 256 * 256) + 0);


    wsPIC(xx, yy - 1, (255 <<24) + (255 * 256 * 256) + (0 * 256 * 256) + 0);
  }
}

vector<ModelTriangle> vertex_normals(vector<ModelTriangle> tris) {
  int ii=-1;
  for (auto i:tris) {
    ii += 1;
    ModelTriangle t = i;
    vector<vec3> normals;
    for (int v = 0; v < t.vert.size(); v++) {
      vec3 vertex = t.normal;
      int count = 1;
      int jj = -1;
      for (auto j:tris) {
        jj += 1;
        ModelTriangle t_ = j;
        for (int u = 0; u < t_.vert.size(); u++)
          if (ii != jj && t.vert[v].x == t_.vert[u].x &&
              t.vert[v].y == t_.vert[u].y &&
              t.vert[v].z == t_.vert[u].z)
            if (acosl(dot(normalize(t.normal), normalize(t_.normal)) /
                      (length(t.normal) * length(t_.normal))) < PI / 4) {
              vertex = vertex + t_.normal;
              count = count + 1;
            }
      }
      vertex = vertex / double(count);
      tris[ii].normals[v] = normalize(vertex);
    }
  }
  return tris;
}



vector<ModelTriangle> parse_obj(string filename, double scale,
                                unordered_map<string, Colour> cols) {
  std::vector<ModelTriangle> tris;
  std::vector<vec3> vert;
  std::vector<vec3> normals;
  std::vector<TexturePoint> texture_points;
  std::string col;
  std::string texture_name;

  std::ifstream File(filename);
  std::string line;

  while (getline(File, line)) {
    if (line == "") continue;

    vector<string> tokens = split(line, ' ');
    assert(tokens.size()!=0);
    auto first = tokens[0];
    if (first == "vn") {
      vec3 nl(stof(tokens[1]), stof(tokens[2]), stof(tokens[3]));
      normals.pb(nl);
    }
    if (first == "v") {
      vec3 vertex(stof(tokens[1]) * scale, stof(tokens[2]) * scale,
                  stof(tokens[3]) * scale);
      vert.pb(vertex);
    }
    if (first == "vt") {
      texture_points.pb(TexturePoint(stof(tokens[1]), stof(tokens[2])));
    }
    if (first == "f") {
      std::vector<string> l1 = split(tokens[1], '/');
      std::vector<string> l2 = split(tokens[2], '/');
      std::vector<string> l3 = split(tokens[3], '/');
      ModelTriangle tri(vert[stoi(*l1.begin()) - 1],
                             vert[stoi(*l2.begin()) - 1],
                             vert[stoi(*l3.begin()) - 1], cols[col]);
      if (normals.size()) {
        tri.normals[0] = normals[stoi(l1[2]) - 1];
        tri.normals[1] = normals[stoi(l2[2]) - 1];
        tri.normals[2] = normals[stoi(l3[2]) - 1];
      }
      tri.normal =
          cross(vec3(tri.vert[1] - (*tri.vert.begin()) ),
                vec3(tri.vert[2] - (*tri.vert.begin()) ));

      if (texture_points.size() && l1[1] != "") {
        tri.texturePoints[0] = texture_points[stoi(l1[1]) - 1];
        tri.texturePoints[1] = texture_points[stoi(l2[1]) - 1];
        tri.texturePoints[2] = texture_points[stoi(l3[1]) - 1];
      }
      tris.pb(tri);
    }
    if (first == "usemtl") {
      col = tokens[1];
    }
  }
  if (!normals.size()) {
    puts("there are no vertex normals in my obj");
    tris = vertex_normals(tris);
  } else
    puts("there WERE vertex normals in my obj");
  File.close();
  return tris;
}

unordered_map<string, Colour> parse_mtl(string filename) {
  unordered_map<string, Colour> cols;
  cols.max_load_factor(0.25);
  string col_name;

  ifstream File(filename);
  string line;

  while (getline(File, line)) {
    if (line == "") continue;

    vector<string> tokens = split(line, ' ');
    auto first = *tokens.begin();
    if (first == "newmtl")
      col_name = tokens[1];
    if (first == "Kd") {
      Colour col(int(stof(tokens[1]) * 255), int(stof(tokens[2]) * 255),
                    int(stof(tokens[3]) * 255));
      cols[col_name] = col;
    }
    if (first == "map_Kd") {
      Colour col = cols[col_name];
      col.name = tokens[1];
      cols[col_name] = col;
    }
  }
  File.close();
  return cols;
}
const auto _0_1_0 = vec3(0.0, 1.0, 0.0);
const auto _1_0_0 = vec3(1.0, 0.0, 0.0);
const auto _0_0_1 = vec3(0.0, 0.0, 1.0);
void look_at() {
  vec3 f = normalize(cam - o);
  vec3 r = normalize(cross(_0_1_0, f));
  vec3 u = normalize(cross(f, r));

  vector<vec3> v = {r,u,f};
  for (int i=0;i<v.size();i++)
    cam_orientation[i] = v[i];
}


void reset_camera() {
  cam = vec3(0.0, 0.0, 4.0);
  light = vec3(0.0, 1.0, 2.0);
  cam_orientation =
      mat3(_1_0_0, _0_1_0, _0_0_1);
}

mat3 rotation_y(double t) {
  return mat3(vec3(cosl(t), 0.0, sinl(t)), _0_1_0,
              vec3(-sinl(t), 0.0, cosl(t)));
}
mat3 rotation_x(double t) {
  return mat3(_1_0_0, vec3(0.0, cosl(t), -sinl(t)),
              vec3(0.0, sinl(t), cosl(t)));
}
mat3 rotation_z(double t) {
  return mat3(vec3(cosl(t), -sinl(t), 0.0), vec3(sinl(t), cosl(t), 0.0),
              _0_0_1);
}

void orbit(bool orb) {
  if (orb) {
    cam = cam * rotation_y(-PI / 180.0);
    cam = cam * rotation_y(-PI / 180.0);
    look_at();
  }
}

function<void(vector<ModelTriangle>, DrawingWindow &)> drawing = draw_raytrace;

void handleEvent(SDL_Event event, DrawingWindow &window) {
  if (event.type == SDL_KEYDOWN) {
    auto kk = event.key.keysym.sym;
    if (kk == SDLK_PAGEDOWN)
      cam.y -= 0.1;
    if(kk == SDLK_PAGEUP)
      cam.y += 0.1;
    if(kk == SDLK_w)
      cam.z -= 0.1;  // cout << "[" << cam.x << "," << cam.y << "," << cam.z <<
                     // "]" << endl;}
    if(kk == SDLK_a)
      cam.x -= 0.1;
    if(kk == SDLK_s)
      cam.z += 0.1;
    if(kk == SDLK_d)
      cam.x += 0.1;
    if(kk == SDLK_q)
      cam = cam * rotation_y(-PI / 180);
    if(kk == SDLK_e)
      cam = cam * rotation_y(PI / 180);
    if(kk == SDLK_EQUALS)
      cam = cam * rotation_x(-PI / 180);
    if(kk == SDLK_MINUS)
      cam = cam * rotation_x(PI / 180);
    if(kk == SDLK_LEFT)
      cam_orientation = cam_orientation * rotation_y(-PI / 180);
    if(kk == SDLK_RIGHT)
      cam_orientation = cam_orientation * rotation_y(PI / 180);
    if(kk == SDLK_UP)
      cam_orientation = cam_orientation * rotation_x(-PI / 180);
    if(kk == SDLK_DOWN)
      cam_orientation = cam_orientation * rotation_x(PI / 180);
    if(kk == SDLK_z)
      cam = cam * rotation_z(-PI / 180);
    if(kk == SDLK_x)
      cam = cam * rotation_z(PI / 180);
    if(kk == SDLK_o)
      orbiting = (orbiting) ? false : true;
    if(kk == SDLK_l)
      look_at();
    if(kk == SDLK_r)
      reset_camera();
    if(kk == SDLK_1) {
      drawing = draw_raytrace;
      cout << "[drawing]: raytrace" << endl;
    } if(kk == SDLK_2) {
      drawing = draw_rasterise;
      cout << "[drawing]: rasterise" << endl;
    } if(kk == SDLK_3) {
      drawing = draw_wireframe;
      cout << "[drawing]: wireframe" << endl;
    } if(kk == SDLK_4) {
      brightness = get_scale;
      cout << "[lighting]: scale" << endl;
    } if(kk == SDLK_5) {
      brightness = gourad;
      cout << "[lighting]: gourad" << endl;
    } if(kk == SDLK_6) {
      brightness = phong;
      cout << "[lighting]: phong" << endl;
    } if(kk == SDLK_KP_8)
      light.z -= 0.1;
    if(kk == SDLK_KP_2)
      light.z += 0.1;
    if(kk == SDLK_KP_6)
      light.x += 0.1;
    if(kk == SDLK_KP_4)
      light.x -= 0.1;
    if(kk == SDLK_KP_MINUS)
      light.y -= 0.1;
    if(kk == SDLK_KP_PLUS)
      light.y += 0.1;
    if(kk == SDLK_LEFTBRACKET) {
      proximity = (proximity) ? false : true;
      cout << "[proximity]: " << proximity << endl;
    } if(kk == SDLK_RIGHTBRACKET) {
      angle_of = (angle_of) ? false : true;
      cout << "[angle_of]: " << angle_of << endl;
    } if(kk == SDLK_HASH) {
      shadows = (shadows) ? false : true;
      cout << "[shadows]: " << shadows << endl;
    } if(kk == SDLK_QUOTE) {
      specular = (specular) ? false : true;
      cout << "[specular]: " << specular << endl;
    } if(kk == SDLK_p)
      show_light = (show_light) ? false : true;


  } if(event.type == SDL_MOUSEBUTTONDOWN)
    window.savePPM("output.ppm");
}

int main(int argc, char *argv[]) {
  vector<ModelTriangle> t =
      parse_obj("cornell-box.obj", 0.5, parse_mtl("cornell-box.mtl"));
  vector<ModelTriangle> tri_2 =
      parse_obj("sphere.obj", 0.5, parse_mtl("cornell-box.mtl"));
  t.insert(t.end(), tri_2.begin(), tri_2.end());

  DrawingWindow window_grey = DrawingWindow(WIDTH, HEIGHT, false);
  SDL_Event event;
  while (true) {
    if (window_grey.pollForInputEvents(event)) handleEvent(event, window_grey);
    orbit(orbiting);

    draw(window_grey);
    drawing(t, window_grey);
    window_grey.renderFrame();
  }
}

