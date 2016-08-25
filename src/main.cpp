#include <SDL2/SDL.h>

#include <iostream>
#include <utility>
#include <cstddef>
#include <cmath>
#include <functional>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <limits>

namespace
{
    constexpr int WINDOW_WIDTH  = 640;
    constexpr int WINDOW_HEIGHT = 480;

    template<typename T>
    class Matrix;

    template<typename T>
    class Vector
    {
        union
        {
            struct {T x, y, z, w;};
            T data[4]{};
        };

    public:
        constexpr Vector() noexcept {};

        constexpr Vector(T x, T y, T z, T w = static_cast<T>(1)) noexcept : data{x, y, z, w} {}

        constexpr Vector& operator+=(const Vector& vector) noexcept {return *this = *this + vector;}

        friend constexpr Vector operator*(const Matrix<T>& matrix, const Vector& vector) noexcept
        {
            Vector temp;
            for (auto i = 0; i < 4; ++i)
                *(temp.data + i) = matrix(i, 0) * *vector.data + matrix(i, 1) * *(vector.data + 1) +
                                                                 matrix(i, 2) * *(vector.data + 2) +
                                                                 matrix(i, 3) * *(vector.data + 3);
            return temp;
        }

        constexpr Vector operator+(const Vector& vector) const noexcept
        {
            return {x + vector.x, y + vector.y, z + vector.z, w + vector.w};
        }

        constexpr Vector operator-(const Vector& vector) const noexcept
        {
            return {x - vector.x, y - vector.y, z - vector.z, w - vector.w};
        }

        friend constexpr Vector operator*(T value, const Vector& vector) noexcept
        {
            return {value * vector.x, value * vector.y, value * vector.z, value * vector.w};
        }

        constexpr Vector operator/(T value) const noexcept {return {x / value, y / value, z / value, w / value};}

        constexpr T get_x() const noexcept {return x;}
        constexpr T get_y() const noexcept {return y;}
        constexpr T get_z() const noexcept {return z;}
        constexpr T get_w() const noexcept {return w;}
    };

    using vec4f = Vector<float>;

    template<typename T>
    class Matrix
    {
        T data[16]{};

    public:
        constexpr Matrix() = default;

        explicit constexpr Matrix(const T (&data)[16]) noexcept
        {
            for (auto i = 0; i < 16; ++i) *(this->data + i) = *(data + i);
        }

        static constexpr Matrix createScreenspace(T width, T height) noexcept
        {
            const auto half_width = width / 2, half_height = height / 2;

            return Matrix{{half_width,          T{},               T{},       half_width,
                                  T{}, -half_height,               T{},      half_height,
                                  T{},          T{}, static_cast<T>(1),              T{},
                                  T{},          T{},               T{}, static_cast<T>(1)}};
        }

        static constexpr Matrix createPerspective(T fov, T aspect, T znear, T zfar) noexcept
        {
            const T f = 1 / std::tan(fov / 2), A =    (zfar + znear) / (znear - zfar),
                                               B = 2 * zfar * znear  / (znear - zfar);

            return Matrix{{f / aspect, T{}, T{}, T{},
                                  T{},   f, T{}, T{},
                                  T{}, T{},   A,   B,
                                  T{}, T{},  -1,  T{}}};
        }

        static constexpr Matrix createScaling(T x, T y, T z) noexcept
        {
            return Matrix{{x, T{}, T{}, T{},
                           T{}, y, T{}, T{},
                           T{}, T{}, z, T{},
                           T{}, T{}, T{}, static_cast<T>(1)}};
        }

        static constexpr Matrix createRotationX(T angle) noexcept
        {
            const T cos_ = std::cos(angle), sin_ = std::sin(angle);

            return Matrix{{static_cast<T>(1),  T{},   T{},              T{},
                                         T{}, cos_, -sin_,              T{},
                                         T{}, sin_,  cos_,              T{},
                                         T{},  T{},   T{}, static_cast<T>(1)}};
        }

        static constexpr Matrix createRotationY(T angle) noexcept
        {
            const T cos_ = std::cos(angle), sin_ = std::sin(angle);

            return Matrix{{ cos_,               T{}, sin_,               T{},
                             T{}, static_cast<T>(1),  T{},               T{},
                           -sin_,               T{}, cos_,               T{},
                             T{},               T{},  T{}, static_cast<T>(1)}};
        }

        static constexpr Matrix createTranslation(T x, T y, T z) noexcept
        {
            return Matrix{{static_cast<T>(1),               T{},               T{},                x,
                                         T{}, static_cast<T>(1),               T{},                y,
                                         T{},               T{}, static_cast<T>(1),                z,
                                         T{},               T{},               T{}, static_cast<T>(1)}};
        }

        constexpr Matrix& operator*=(const Matrix& matrix) noexcept {return *this = *this * matrix;}

        constexpr Matrix operator*(const Matrix& matrix) const noexcept
        {
            Matrix temp;
            for (auto i = 0; i < 4; ++i)
                for (auto j = 0; j < 4; ++j)
                    for (auto k = 0; k < 4; ++k)
                        *(temp.data + i * 4 + j) += *(data + i * 4 + k) * *(matrix.data + k * 4 + j);

            return temp;
        }

        constexpr T operator()(unsigned row, unsigned col) const noexcept {return *(data + 4 * row + col);}
    };

    using mat4f = Matrix<float>;

    struct Vertex
    {
        vec4f position;
        vec4f color;
    };

    struct Edge
    {
        Vertex vert1;
        Vertex vert2;
    };

    struct Mesh
    {
        std::vector<Vertex> vertices;
    };

    float zbuffer[WINDOW_WIDTH * WINDOW_HEIGHT];

    void clear_zbuffer() noexcept
    {
        for (size_t i = 0; i < WINDOW_WIDTH * WINDOW_HEIGHT; ++i)
            *(zbuffer + i) = std::numeric_limits<float>::max();
    }

    Mesh load_obj(const std::string& filepath)
    {
        std::ifstream stream{filepath, std::ios::in};
        if (!stream)
            throw std::runtime_error{"obj loading error"};

        std::vector<vec4f> vertex_positions, normals;
        std::vector<size_t> vertex_indices, normal_indices;
        std::string line;

        while (std::getline(stream, line))
        {
            std::stringstream line_stream{line};
            std::string c;
            line_stream >> c;

            if (c == "v")
            {
                float x, y, z;
                line_stream >> x >> y >> z;
                vertex_positions.emplace_back(x, y, z);
            }
            else if (c == "vn")
            {
                float x, y, z;
                line_stream >> x >> y >> z;
                normals.emplace_back(x, y, z);
            }
            else if (c == "f")
            {
                for (auto i = 0; i < 3; ++i)
                {
                    size_t vertex_index;
                    line_stream >> vertex_index;
                    line_stream.ignore(2);

                    size_t normal_index;
                    line_stream >> normal_index;

                    vertex_indices.push_back(vertex_index);
                    normal_indices.push_back(normal_index);
                }
            }
        }

        std::vector<Vertex> vertices;

        for (auto i : vertex_indices)
            vertices.push_back({vertex_positions[i - 1], normals[normal_indices[i - 1]]});

        return {vertices};
    }

    inline void set_pixel(SDL_Surface* surface, int x, int y, Uint8 r, Uint8 g, Uint8 b) noexcept
    {
        *(static_cast<Uint32*>(surface->pixels) + surface->w * y + x) = ::SDL_MapRGB(surface->format, r, g, b);
    }

    template<typename ComparisonFunc>
    void scanline(SDL_Surface* surface, const Edge& edge1, const Edge& edge2, ComparisonFunc compare) noexcept
    {
        const vec4f& pos1 = edge1.vert1.position, pos2 = edge1.vert2.position,
                     pos3 = edge2.vert1.position, pos4 = edge2.vert2.position;

        const auto pos21 = pos2 - pos1, pos43 = pos4 - pos3;

        if (!pos21.get_y())
            return;

        auto ldx = pos21.get_x() / pos21.get_y();
        auto rdx = pos43.get_x() / pos43.get_y();

        auto ldz = pos21.get_z() / pos21.get_y();
        auto rdz = pos43.get_z() / pos43.get_y();

        auto lx = pos1.get_x();
        auto rx = pos3.get_x() + (pos1.get_y() - pos3.get_y()) * rdx;

        auto lz = pos1.get_z();
        auto rz = pos3.get_z() + (pos1.get_y() - pos3.get_y()) * rdz;

        const vec4f& color1 = edge1.vert1.color, color2 = edge1.vert2.color,
                     color3 = edge2.vert1.color, color4 = edge2.vert2.color;

        const auto color21 = color2 - color1, color43 = color4 - color3;
        
        auto ldc = color21 / pos21.get_y();
        auto rdc = color43 / pos43.get_y();

        auto lc = color1;
        auto rc = color3 + (pos1.get_y() - pos3.get_y()) * rdc;

        if (compare(ldx, rdx))
        {
            std::swap(ldx, rdx); std::swap(rx, lx);
            std::swap(ldc, rdc); std::swap(rc, lc);
        }

        const int min_y = std::ceil(pos1.get_y()), max_y = std::ceil(pos2.get_y());
        if (min_y < 0 || max_y >= surface->h)
            return;

        for (int y = min_y; y < max_y; ++y)
        {
            const int min_x = std::ceil(lx), max_x = std::ceil(rx);
            if (min_x < 0 || max_x >= surface->w)
                return;

            const auto pixel_delta_color = (rc - lc) / (rx - lx);
            auto pixel_color = lc;

            for (int x = min_x; x < max_x; ++x)
            {
                const auto zvalue = lz + (x - min_x) * (rz - lz) / (rx - lx);

                if (zvalue < *(zbuffer + WINDOW_WIDTH * y + x))
                {
                    *(zbuffer + WINDOW_WIDTH * y + x) = zvalue;
                    ::set_pixel(surface, x, y, pixel_color.get_x() * 255,
                                               pixel_color.get_y() * 255,
                                               pixel_color.get_z() * 255);
                }
                pixel_color += pixel_delta_color;
            }

            lx += ldx;
            rx += rdx;
            lz += ldz;
            rz += rdz;
            lc += ldc;
            rc += rdc;
        }
    }

    void rasterize_triangle(SDL_Surface* surface, Vertex vert1, Vertex vert2, Vertex vert3) noexcept
    {
        if (vert1.position.get_y() > vert2.position.get_y())
            std::swap(vert1, vert2);
        if (vert1.position.get_y() > vert3.position.get_y())
            std::swap(vert1, vert3);
        if (vert2.position.get_y() > vert3.position.get_y())
            std::swap(vert2, vert3);

        ::scanline(surface, {vert1, vert2}, {vert1, vert3}, std::greater<>{});
        ::scanline(surface, {vert2, vert3}, {vert1, vert3}, std::less   <>{});
    }

    inline constexpr bool check_virtual_position(const vec4f& position) noexcept
    {
        return position.get_x() < -position.get_w() || position.get_x() > position.get_w() ||
               position.get_y() < -position.get_w() || position.get_y() > position.get_w() ||
               position.get_z() < -position.get_w() || position.get_z() > position.get_w();
    }

    void render_triangle(SDL_Surface* surface, const Vertex& vert1, const Vertex& vert2,
                                                                    const Vertex& vert3, const mat4f& matrix) noexcept
    {
        const vec4f virtual_pos1 = matrix * vert1.position, virtual_pos2 = matrix * vert2.position,
                    virtual_pos3 = matrix * vert3.position;

        if (::check_virtual_position(virtual_pos1) ||
            ::check_virtual_position(virtual_pos2) || ::check_virtual_position(virtual_pos3)) return;

        const auto screenspace_matrix = mat4f::createScreenspace(surface->w, surface->h);

        const auto screenspace_pos1 = screenspace_matrix * virtual_pos1,
                   screenspace_pos2 = screenspace_matrix * virtual_pos2,
                   screenspace_pos3 = screenspace_matrix * virtual_pos3;

        ::rasterize_triangle(surface, {screenspace_pos1 / screenspace_pos1.get_w(), vert1.color},
                                      {screenspace_pos2 / screenspace_pos2.get_w(), vert2.color},
                                      {screenspace_pos3 / screenspace_pos3.get_w(), vert3.color});
    }

    void render_mesh(SDL_Surface* surface, const Mesh& mesh, const mat4f& matrix) noexcept
    {
        for (decltype(mesh.vertices)::size_type i = 0; i < mesh.vertices.size(); i += 3)
            ::render_triangle(surface, mesh.vertices[i], mesh.vertices[i + 1], mesh.vertices[i + 2], matrix);
    }
}

int main()
{
    if (::SDL_Init(SDL_INIT_VIDEO) >= 0)
    {
        SDL_Window* const window = ::SDL_CreateWindow("Software Renderer",
                SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, WINDOW_WIDTH, WINDOW_HEIGHT, SDL_WINDOW_SHOWN);
        if (window)
        {
            SDL_Surface* const surface = ::SDL_GetWindowSurface(window);
            if (surface)
            {
                try
                {
                    const Mesh mesh{::load_obj("res/bunny.obj")};

                    const auto perspective_matrix = mat4f::createPerspective(M_PI / 3,
                            static_cast<float>(WINDOW_WIDTH) / WINDOW_HEIGHT, 1e-1F, 1e3F);
                    auto matrix1 = perspective_matrix * mat4f::createTranslation(0.0F, -0.12F, -0.25F);
                    auto matrix2 = perspective_matrix * mat4f::createTranslation(2.0F, 1.0F, -5.0F) *
                            mat4f::createScaling(0.4F, 0.4F, 0.4F);
                    const auto rotation_matrix = mat4f::createRotationY(0.05F) * mat4f::createRotationX(-0.01F);

                    SDL_Event event;

                    bool running = true;
                    while (running)
                    {
                        while (::SDL_PollEvent(&event))
                        {
                            if (event.type == SDL_QUIT)
                                running = false;
                        }

                        ::clear_zbuffer();

                        ::SDL_FillRect(surface, nullptr, 0);
                        ::SDL_LockSurface(surface);
                        ::render_mesh(surface, mesh, matrix1);
                        ::render_triangle(surface, {{ 0.0F, -1.0F, 0.0F}, {1.0F, 0.0F, 0.0F}},
                                                   {{-1.0F,  1.0F, 0.0F}, {0.0F, 1.0F, 0.0F}},
                                                   {{ 1.0F,  1.0F, 0.0F}, {0.0F, 0.0F, 1.0F}}, matrix2);
                        ::SDL_UnlockSurface(surface);

                        matrix1 *= rotation_matrix;
                        matrix2 *= rotation_matrix;

                        ::SDL_UpdateWindowSurface(window);
                        ::SDL_Delay(20);
                    }
                }
                catch (const std::exception& ex)
                {
                    std::cerr << ex.what() << std::endl;
                }
            }
            else
                std::cerr << "SDL_Surface getting error: " << ::SDL_GetError() << std::endl;

            ::SDL_DestroyWindow(window);
        }
        else
            std::cerr << "SDL_Window creation error: " << ::SDL_GetError() << std::endl;

        ::SDL_Quit();
    }
    else
        std::cerr << "SDL2 initialization error: " << ::SDL_GetError() << std::endl;

    return 0;
}
