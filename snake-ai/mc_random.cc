#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <list>
#include <stack>
#include <deque>
#include <queue>
#include <vector>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "jsoncpp/json.h"

using std::cin;
using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::queue;
using std::deque;
using std::stack;
using std::list;
using std::string;

const int MAXN = 27;
const int INF = 0x3F3F3F3F;

struct Point {
    int x, y;
    Point(int a, int b): x(a), y(b) {}
    Point& operator+=(const Point& rhs) { x += rhs.x; y += rhs.y; return *this; }
    Point& operator-=(const Point& rhs) { x -= rhs.x; y -= rhs.y; return *this; }
    friend inline Point operator+(Point lhs, const Point& rhs) { return lhs += rhs; }
    friend inline Point operator-(Point lhs, const Point& rhs) { return lhs -= rhs; }
    friend inline bool operator==(const Point& lhs, const Point& rhs) { return lhs.x == rhs.x && lhs.y == rhs.y; }
};
const Point XY[4] = {Point(-1, 0), Point(0, 1), Point(1, 0), Point(0, -1)};

int n, m;
int cur_round;
char map[MAXN][MAXN];
char fmap[MAXN][MAXN];
std::ostringstream debug;
clock_t clock_start = clock();

inline double get_time() {
    return static_cast<double>(clock()-clock_start)/CLOCKS_PER_SEC;
}

int dist(const Point& p1, const Point& p2) {
    return abs(p1.x-p2.x) + abs(p1.y-p2.y);
}


inline bool will_grow(int round) {
    if (round <= 9) return true;
    return (round-9) % 3 == 0;
}

struct Snake {
    int id;
    char (&local_map)[MAXN][MAXN];
    deque<Point> body;
    
    Snake(): id(0), local_map(::map), body() {}
    
    Snake(const int i, char (&m)[MAXN][MAXN], const deque<Point> &b):
    id(i), local_map(m), body(b) {}
    
    Snake fork(char (&m)[MAXN][MAXN]) const {
        return Snake(id, m, body);
    }
    
    
    bool in_body(int x, int y) const {
        for (deque<Point>::const_iterator iter = body.begin(); iter != body.end(); ++iter)
            if (x == iter->x && y == iter->y)
                return true;
        return false;
    }
    
    void insert_front(const Point& p) {
        body.push_front(p);
        local_map[p.x][p.y] = '0' + id;
    }
    
    void insert_back(const Point& p) {
        body.push_back(p);
        local_map[p.x][p.y] = '0' + id;
    }
    
    void delete_front() {
        const Point p(body.front());
        local_map[p.x][p.y] = '.';
        body.pop_front();
    }
    
    void delete_back() {
        const Point p(body.back());
        local_map[p.x][p.y] = '.';
        body.pop_back();
    }
    
    void print() const {
        for (deque<Point>::const_iterator iter = body.begin(); iter != body.end(); ++iter)
            printf("(%d, %d) ", iter->x, iter->y);
        puts("");
    }
    
    bool valid_direction(int dir) const {
        const Point p = body.front() + XY[dir];
        if (p.x > n || p.y > m || p.x < 1 || p.y < 1) return false;
        if (local_map[p.x][p.y] != '.') return false;
        return true;
    }
    
    int dist(const Point& p) const {
        int d = 0x7FFFFFFF;
        for (deque<Point>::const_iterator iter = body.begin(); iter != body.end(); ++iter)
            d = std::min(d, ::dist(*iter, p));
        return d;
    }
    
    int body_leave_round(const Point& p) const {
        int round = cur_round+1;
        for (deque<Point>::const_reverse_iterator iter = body.rbegin(); iter != body.rend(); ++iter) {
            if (will_grow(round)) ++round;
            if (*iter == p) return round;
            ++round;
        }
        return INF;
    }
} snake[2];

void print_map(const char (&map)[MAXN][MAXN] = ::map, std::ostream& out = cerr) {
    out << "===========print_map start=======" << endl;
    for (int j = 1; j <= m; ++j) {
        for (int i = 1; i <= n; ++i)
            if (Point(i, j) == snake[0].body.front()) out << '@';
            else if (Point(i, j) == snake[1].body.front()) out << '%';
            else out << map[i][j];
        out << endl;
    }
    out << "===========print_map end=========" << endl;
}

void print_array(const int (&a)[MAXN][MAXN], std::ostream& out = cout) {
    out << "===========print_array start====" << endl;
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j)
            out << a[i][j] << ' ';
        out << endl;
    }
    out << "===========print_array end======" << endl;
}



void init() {
    string str, tmp;
    while (getline(cin, tmp)) str += tmp;
    
    Json::Reader reader;
    Json::Value input;
    reader.parse(str, input);
    n = input["requests"][static_cast<Json::Value::UInt>(0)]["height"].asInt();
    m = input["requests"][static_cast<Json::Value::UInt>(0)]["width"].asInt();
    
    for (int i = 1; i <= n; ++i)
        for (int j = 1; j <= m; ++j)
            map[i][j] = '.';
    
    const int x = input["requests"][static_cast<Json::Value::UInt>(0)]["x"].asInt();
    snake[0].id = 0, snake[1].id = 1;
    if (x==1) {
        snake[0].insert_front(Point(1,1));
        snake[1].insert_front(Point(n,m));
    } else {
        snake[1].insert_front(Point(1,1));
        snake[0].insert_front(Point(n,m));
    }
    
    const int cnt_obstacle = input["requests"][static_cast<Json::Value::UInt>(0)]["obstacle"].size();
    for (int i = 0; i < cnt_obstacle; ++i) {
        const int x = input["requests"][static_cast<Json::Value::UInt>(0)]["obstacle"][static_cast<Json::Value::UInt>(i)]["x"].asInt();
        const int y = input["requests"][static_cast<Json::Value::UInt>(0)]["obstacle"][static_cast<Json::Value::UInt>(i)]["y"].asInt();
        map[x][y] = '#';
    }
    
    const int cnt_round = input["responses"].size();
    for (int i = 0; i < cnt_round; ++i) {
        const int dir0 = input["responses"][i]["direction"].asInt();
        const int dir1 = input["requests"][i+1]["direction"].asInt();
        
        if (!will_grow(i)) {
            snake[0].delete_back();
            snake[1].delete_back();
        }
        snake[0].insert_front(snake[0].body.front()+XY[dir0]);
        snake[1].insert_front(snake[1].body.front()+XY[dir1]);
    }
    
    cur_round = cnt_round;
    if (!will_grow(cnt_round)) {
        snake[0].delete_back();
        snake[1].delete_back();
    }
    
    srand(static_cast<unsigned>(time(0)) + cnt_round);
}

Json::Value ret;
void write_result(const int res) {
    debug << "time=" << get_time() << endl;
    ret["response"]["direction"] = res;
    ret["response"]["debug"] = debug.str();
    Json::FastWriter writer;
    cout << writer.write(ret) << endl;
    exit(0);
}

const int DETERMINED_SUCCESS = 10;
const int DETERMINED_FAILURE = -10;
const int UNDETERMINED = 0;
inline void recover_tail(const int round, Snake& self, Snake& oppo, const Point& tail0, const Point& tail1) {
    if (!will_grow(round)) {
        self.insert_back(tail0);
        oppo.insert_back(tail1);
    }
}
int check_determined(const int round, const int step, Snake& self, Snake& oppo, int sol[] = 0) {
    if (step == 14) return 0;
    const Point tail0(self.body.back());
    const Point tail1(oppo.body.back());
    if (!will_grow(round)) {
        self.delete_back();
        oppo.delete_back();
    }
    int p[4] = {0,0,0,0};
    for (int d = 0; d < 4; ++d) {
        if (!self.valid_direction(d)) continue;
        self.insert_front(self.body.front()+XY[d]);
        const int ret = check_determined(round+1, step+1, oppo, self);
        self.delete_front();
        
        if (ret == DETERMINED_FAILURE) {
            recover_tail(round, self, oppo, tail0, tail1);
            return DETERMINED_SUCCESS + d;
        }
        if (ret >= DETERMINED_SUCCESS) p[d] = 1;
        else p[d] = 100;
    }
    recover_tail(round, self, oppo, tail0, tail1);
    if (!p[0] && !p[1] && !p[2] && !p[3]) return DETERMINED_FAILURE;
    if (sol) memcpy(sol, p, sizeof(p));
    return UNDETERMINED;
}

inline int choose_direction(const int dtm, const int (&p)[4]) {
    if (dtm == DETERMINED_FAILURE) return -1;
    if (dtm == DETERMINED_SUCCESS) return dtm-DETERMINED_SUCCESS;
    
    // random choose
    const int s0 = p[0];
    const int s1 = s0 + p[1];
    const int s2 = s1 + p[2];
    const int s3 = s2 + p[3];
    const int r = rand() % s3;
    if (s0 > r) return 0;
    if (s1 > r) return 1;
    if (s2 > r) return 2;
    return 3;
}

void decide() {
    int score[4] = {0, 0, 0, 0};
    int cnt_random = 0;
    
    bool any = false;
    for (int rd1 = 0; rd1 < 4; ++rd1)
        any = any || snake[1].valid_direction(rd1);
    if (!any) {
        for (int rd0 = 0; rd0 < 4; ++rd0)
            if (snake[0].valid_direction(rd0))
                write_result(rd0);
        write_result(-1);
    }
    
    while (get_time() < 0.8) {
        for (int rd0 = 0; rd0 < 4; ++rd0) {
            if (!snake[0].valid_direction(rd0)) continue;
            for (int rd1 = 0; rd1 < 4; ++rd1) {
                if (!snake[1].valid_direction(rd1)) continue;
                
                memcpy(fmap, map, sizeof(map));
                Snake s0(snake[0].fork(fmap));
                Snake s1(snake[1].fork(fmap));
                s0.insert_front(s0.body.front()+XY[rd0]);
                s1.insert_front(s1.body.front()+XY[rd1]);
                
                for (int round = cur_round + 1; ; ++round) {
                    
                    int p0[4], p1[4];
                    const int dtm0 = check_determined(round, 0, s0, s1, p0);
                    const int dtm1 = check_determined(round, 0, s1, s0, p1);
                    const int d0 = choose_direction(dtm0, p0);
                    const int d1 = choose_direction(dtm1, p1);
                    
                    if (d0 != -1 && d1 != -1) {
                        s0.insert_front(s0.body.front()+XY[d0]);
                        s1.insert_front(s1.body.front()+XY[d1]);
                        continue;
                    }
                    
                    ++cnt_random;
                    if (d0 == -1 && d1 == -1)
                        score[rd0] += 1;
                    else if (d1 == -1)
                        score[rd0] += 2;
                    break;
                }
            }
        }
    }
    
    int d = 0;
    debug << "cnt_random = " << cnt_random << ", score = [" << std::fixed << std::setprecision(3);
    for (int i = 0; i < 4; ++i) {
        if (score[i] > score[d])
            d = i;
        debug << score[i] / static_cast<double>(cnt_random) << ',';
    }
    debug << std::resetiosflags( std::ios::fixed | std::ios::showpoint ) << ']' << endl;
    
    if (score[d] == 0) {
        std::vector<int> v;
        for (int i = 0; i < 4; ++i)
            if (snake[0].valid_direction(i))
                v.push_back(i);
        if (!v.empty())
            d = v[rand()%v.size()];
    }
    write_result(d);
}

int main() {
    //freopen("input.json", "r", stdin);
    init();
    decide();
}
