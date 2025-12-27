#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <memory>
#include <random>
#include <fstream>
#include <omp.h> // 需要开启 OpenMP 支持

//Params
const double infinity = std::numeric_limits<double>::infinity();
const double pi = 3.1415926535897932385;

// --- 1. Utility Functions ---
inline double degrees_to_radians(double degrees) {
    return degrees * pi / 180.0;
}

inline double random_double() {
    static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    static std::mt19937 generator;
    return distribution(generator);
}

inline double random_double(double min, double max) {
    return min + (max - min) * random_double();
}

inline double clamp(double x, double min, double max) {
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

// --- 2. Vec3 Class (Math Core) ---
class vec3 {
public:
    double e[3];
    vec3() : e{0,0,0} {}
    vec3(double e0, double e1, double e2) : e{e0, e1, e2} {}

    double x() const { return e[0]; }
    double y() const { return e[1]; }
    double z() const { return e[2]; }

    vec3 operator-() const { return vec3(-e[0], -e[1], -e[2]); }
    vec3& operator+=(const vec3 &v) {
        e[0] += v.e[0]; e[1] += v.e[1]; e[2] += v.e[2];
        return *this;
    }
    vec3& operator*=(const double t) {
        e[0] *= t; e[1] *= t; e[2] *= t;
        return *this;
    }
    double length_squared() const { return e[0]*e[0] + e[1]*e[1] + e[2]*e[2]; }
    double length() const { return std::sqrt(length_squared()); }
    
    static vec3 random() { return vec3(random_double(), random_double(), random_double()); }
    static vec3 random(double min, double max) { return vec3(random_double(min,max), random_double(min,max), random_double(min,max)); }
};

using point3 = vec3;
using color = vec3;

inline std::ostream& operator<<(std::ostream &out, const vec3 &v) {
    return out << v.e[0] << ' ' << v.e[1] << ' ' << v.e[2];
}
inline vec3 operator+(const vec3 &u, const vec3 &v) { return vec3(u.e[0]+v.e[0], u.e[1]+v.e[1], u.e[2]+v.e[2]); }
inline vec3 operator-(const vec3 &u, const vec3 &v) { return vec3(u.e[0]-v.e[0], u.e[1]-v.e[1], u.e[2]-v.e[2]); }
inline vec3 operator*(const vec3 &u, const vec3 &v) { return vec3(u.e[0]*v.e[0], u.e[1]*v.e[1], u.e[2]*v.e[2]); }
inline vec3 operator*(double t, const vec3 &v) { return vec3(t*v.e[0], t*v.e[1], t*v.e[2]); }
inline vec3 operator*(const vec3 &v, double t) { return t * v; }
inline vec3 operator/(vec3 v, double t) { return (1/t) * v; }
inline double dot(const vec3 &u, const vec3 &v) { return u.e[0] * v.e[0] + u.e[1] * v.e[1] + u.e[2] * v.e[2]; }
inline vec3 cross(const vec3 &u, const vec3 &v) {
    return vec3(u.e[1] * v.e[2] - u.e[2] * v.e[1],
                u.e[2] * v.e[0] - u.e[0] * v.e[2],
                u.e[0] * v.e[1] - u.e[1] * v.e[0]);
}
inline vec3 unit_vector(vec3 v) { return v / v.length(); }

vec3 random_in_unit_sphere() {
    while (true) {
        auto p = vec3::random(-1,1);
        if (p.length_squared() >= 1) continue;
        return p;
    }
}
vec3 random_unit_vector() { return unit_vector(random_in_unit_sphere()); }

// --- 3. Ray Class ---
class ray {
public:
    point3 orig;
    vec3 dir;
    ray() {}
    ray(const point3& origin, const vec3& direction) : orig(origin), dir(direction) {}
    point3 origin() const { return orig; }
    vec3 direction() const { return dir; }
    point3 at(double t) const { return orig + t*dir; }
};

// --- 4. Hitable Abstraction ---
class material;

struct hit_record {
    point3 p;
    vec3 normal;
    std::shared_ptr<material> mat_ptr;
    double t;
    bool front_face;

    inline void set_face_normal(const ray& r, const vec3& outward_normal) {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }
};

class hitable {
public:
    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
};

// --- 5. Materials ---
class material {
public:
    virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered) const = 0;
};

class lambertian : public material {
public:
    color albedo;
    lambertian(const color& a) : albedo(a) {}
    virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered) const override {
        auto scatter_direction = rec.normal + random_unit_vector();
        if (scatter_direction.length_squared() < 1e-8) scatter_direction = rec.normal; // Catch degenerate scatter
        scattered = ray(rec.p, scatter_direction);
        attenuation = albedo;
        return true;
    }
};

class metal : public material {
public:
    color albedo;
    double fuzz;
    metal(const color& a, double f) : albedo(a), fuzz(f < 1 ? f : 1) {}
    virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered) const override {
        vec3 reflected = unit_vector(r_in.direction()) - 2*dot(unit_vector(r_in.direction()), rec.normal)*rec.normal;
        scattered = ray(rec.p, reflected + fuzz*random_in_unit_sphere());
        attenuation = albedo;
        return (dot(scattered.direction(), rec.normal) > 0);
    }
};

class dielectric : public material {
public:
    double ir; // Index of Refraction
    dielectric(double index_of_refraction) : ir(index_of_refraction) {}

    virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered) const override {
        attenuation = color(1.0, 1.0, 1.0);
        double refraction_ratio = rec.front_face ? (1.0/ir) : ir;

        vec3 unit_direction = unit_vector(r_in.direction());
        double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
        double sin_theta = sqrt(1.0 - cos_theta*cos_theta);

        bool cannot_refract = refraction_ratio * sin_theta > 1.0;
        vec3 direction;

        if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double())
            direction = unit_vector(r_in.direction()) - 2*dot(unit_vector(r_in.direction()), rec.normal)*rec.normal; // Reflect
        else {
            vec3 r_out_perp =  refraction_ratio * (unit_direction + cos_theta*rec.normal);
            vec3 r_out_parallel = -sqrt(fabs(1.0 - r_out_perp.length_squared())) * rec.normal;
            direction = r_out_perp + r_out_parallel; // Refract
        }

        scattered = ray(rec.p, direction);
        return true;
    }
private:
    static double reflectance(double cosine, double ref_idx) {
        // Schlick's approximation
        auto r0 = (1-ref_idx) / (1+ref_idx);
        r0 = r0*r0;
        return r0 + (1-r0)*pow((1 - cosine),5);
    }
};

// --- 6. Objects (Sphere & World) ---
class sphere : public hitable {
public:
    point3 center;
    double radius;
    std::shared_ptr<material> mat_ptr;

    sphere() {}
    sphere(point3 cen, double r, std::shared_ptr<material> m) : center(cen), radius(r), mat_ptr(m) {}

    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
        vec3 oc = r.origin() - center;
        auto a = r.direction().length_squared();
        auto half_b = dot(oc, r.direction());
        auto c = oc.length_squared() - radius*radius;
        auto discriminant = half_b*half_b - a*c;
        if (discriminant < 0) return false;
        auto sqrtd = sqrt(discriminant);

        auto root = (-half_b - sqrtd) / a;
        if (root < t_min || root > t_max) {
            root = (-half_b + sqrtd) / a;
            if (root < t_min || root > t_max) return false;
        }

        rec.t = root;
        rec.p = r.at(rec.t);
        vec3 outward_normal = (rec.p - center) / radius;
        rec.set_face_normal(r, outward_normal);
        rec.mat_ptr = mat_ptr;
        return true;
    }
};

class hitable_list : public hitable {
public:
    std::vector<std::shared_ptr<hitable>> objects;
    hitable_list() {}
    void add(std::shared_ptr<hitable> object) { objects.push_back(object); }

    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override {
        hit_record temp_rec;
        bool hit_anything = false;
        auto closest_so_far = t_max;

        for (const auto& object : objects) {
            if (object->hit(r, t_min, closest_so_far, temp_rec)) {
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec;
            }
        }
        return hit_anything;
    }
};

// --- 7. Camera ---
class camera {
public:
    point3 origin;
    point3 lower_left_corner;
    vec3 horizontal;
    vec3 vertical;

    camera(point3 lookfrom, point3 lookat, vec3 vup, double vfov, double aspect_ratio) {
        auto theta = degrees_to_radians(vfov);
        auto h = tan(theta/2);
        auto viewport_height = 2.0 * h;
        auto viewport_width = aspect_ratio * viewport_height;

        auto w = unit_vector(lookfrom - lookat);
        auto u = unit_vector(cross(vup, w));
        auto v = cross(w, u);

        origin = lookfrom;
        horizontal = viewport_width * u;
        vertical = viewport_height * v;
        lower_left_corner = origin - horizontal/2 - vertical/2 - w;
    }

    ray get_ray(double s, double t) const {
        return ray(origin, lower_left_corner + s*horizontal + t*vertical - origin);
    }
};

// --- 8. The Renderer ---
color ray_color(const ray& r, const hitable& world, int depth) {
    hit_record rec;

    // Bounce limit exceeded, no more light
    if (depth <= 0) return color(0,0,0);

    // 0.001 to ignore hits very near zero (floating point shadow acne)
    if (world.hit(r, 0.001, infinity, rec)) {
        ray scattered;
        color attenuation;
        if (rec.mat_ptr->scatter(r, rec, attenuation, scattered))
            return attenuation * ray_color(scattered, world, depth-1);
        return color(0,0,0);
    }
    
    // Background (Gradient Blue)
    vec3 unit_direction = unit_vector(r.direction());
    auto t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*color(1.0, 1.0, 1.0) + t*color(0.5, 0.7, 1.0);
}

void write_color(std::ostream &out, color pixel_color, int samples_per_pixel) {
    auto r = pixel_color.x();
    auto g = pixel_color.y();
    auto b = pixel_color.z();

    // Divide the color by the number of samples and gamma-correct for gamma=2.0.
    auto scale = 1.0 / samples_per_pixel;
    r = sqrt(scale * r);
    g = sqrt(scale * g);
    b = sqrt(scale * b);

    // Write the translated [0,255] value of each color component.
    out << static_cast<int>(256 * clamp(r, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(g, 0.0, 0.999)) << ' '
        << static_cast<int>(256 * clamp(b, 0.0, 0.999)) << '\n';
}

int main() {
    // Image
    const auto aspect_ratio = 16.0 / 9.0;
    const int image_width = 800; // Increase for higher res
    const int image_height = static_cast<int>(image_width / aspect_ratio);
    const int samples_per_pixel = 50; // Higher = slower but less noise
    const int max_depth = 50;

    // World
    hitable_list world;
    auto material_ground = std::make_shared<lambertian>(color(0.8, 0.8, 0.0));
    auto material_center = std::make_shared<lambertian>(color(0.1, 0.2, 0.5));
    auto material_left   = std::make_shared<dielectric>(1.5);
    auto material_right  = std::make_shared<metal>(color(0.8, 0.6, 0.2), 0.0);

    world.add(std::make_shared<sphere>(point3( 0.0, -100.5, -1.0), 100.0, material_ground));
    world.add(std::make_shared<sphere>(point3( 0.0,    0.0, -1.0),   0.5, material_center));
    world.add(std::make_shared<sphere>(point3(-1.0,    0.0, -1.0),   0.5, material_left));
    world.add(std::make_shared<sphere>(point3(-1.0,    0.0, -1.0), -0.45, material_left)); // Hollow glass
    world.add(std::make_shared<sphere>(point3( 1.0,    0.0, -1.0),   0.5, material_right));

    // Camera
    point3 lookfrom(3,3,2);
    point3 lookat(0,0,-1);
    vec3 vup(0,1,0);
    auto dist_to_focus = (lookfrom-lookat).length();
    auto aperture = 2.0;
    camera cam(lookfrom, lookat, vup, 20, aspect_ratio);

    // Render
    std::cout << "P3\n" << image_width << " " << image_height << "\n255\n";

    // Buffer to store pixels (for parallel processing)
    std::vector<color> buffer(image_width * image_height);

    // Parallel Rendering Loop
    #pragma omp parallel for schedule(dynamic)
    for (int j = 0; j < image_height; ++j) {
        // Optional: Progress indicator (thread unsafe output, just for rough idea)
        // if (j % 50 == 0) std::cerr << "\rScanlines remaining: " << (image_height - j) << ' ' << std::flush;
        
        for (int i = 0; i < image_width; ++i) {
            color pixel_color(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                auto u = (i + random_double()) / (image_width-1);
                auto v = (image_height - 1 - j + random_double()) / (image_height-1); // Flip V for image coords
                ray r = cam.get_ray(u, v);
                pixel_color += ray_color(r, world, max_depth);
            }
            buffer[j * image_width + i] = pixel_color;
        }
    }

    // Output Image
    for (const auto& pixel : buffer) {
        write_color(std::cout, pixel, samples_per_pixel);
    }

    std::cerr << "\nDone.\n";
}