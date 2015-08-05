#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <list>
#include <deque>
#include <queue>
#include <vector>
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
using std::list;
using std::string;

const int MAXN = 27;
const int INF = 0x3F3F3F3F;

const int DIE_STEP = 10;

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
int board_dist[MAXN][MAXN];
int board_from[MAXN][MAXN];
char map[MAXN][MAXN];
int cnt_sol_not_die;
bool follow_sol_not_die;
char fmap[MAXN][MAXN];
deque<int> sol_not_die;
queue<Point> Q;
std::ostringstream debug;

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
    list<Point> body;

    Snake(): id(0), local_map(::map), body() {}

    Snake(const int i, char (&m)[MAXN][MAXN], const list<Point> &b):
        id(i), local_map(m), body(b) {}

    Snake fork(char (&m)[MAXN][MAXN]) const {
        return Snake(id, m, body);
    }


    bool in_body(int x, int y) const {
        for (list<Point>::const_iterator iter = body.begin(); iter != body.end(); ++iter)
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

    void move(int dir, int round) {
        const Point &p(body.front());
        if (!will_grow(round)) delete_back();
        insert_front(Point(p+XY[dir]));
    }

    void print() const {
        for (list<Point>::const_iterator iter = body.begin(); iter != body.end(); ++iter)
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
        for (list<Point>::const_iterator iter = body.begin(); iter != body.end(); ++iter)
            d = std::min(d, ::dist(*iter, p));
        return d;
    }
    
    int body_leave_round(const Point& p) const {
        int round = cur_round+1;
        for (list<Point>::const_reverse_iterator iter = body.rbegin(); iter != body.rend(); ++iter) {
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

void print_snake_map(const Snake& s) {
    for (int j = 1; j <= m; ++j) {
        for (int i = 1; i <= n; ++i)
            if (Point(i, j) == s.body.front()) cerr << '@';
            else cerr << map[i][j];
        cerr << endl;
    }
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
        snake[0].move(dir0, i);
        snake[1].move(dir1, i);
    }
    
    follow_sol_not_die = input["responses"][cnt_round-1]["data"]["follow_sol_not_die"].asBool();
    debug << "follow_sol_not_die=" << follow_sol_not_die << endl;
    debug << "sol_not_die=[";
    int cnt = input["responses"][cnt_round-1]["data"]["sol_not_die"].size();
    for (int i = 0; i < cnt; ++i) {
        sol_not_die.push_back(input["responses"][cnt_round-1]["data"]["sol_not_die"][i].asInt());
        debug << sol_not_die.back() << ',';
    }
    debug << ']' << endl;
    

    cur_round = cnt_round;
    if (!will_grow(cnt_round)) {
        snake[0].delete_back();
        snake[1].delete_back();
    }

    srand(static_cast<unsigned>(time(0)) + cnt_round);
}

Json::Value ret;
void write_result(const int res) {
    ret["response"]["direction"] = res;
    ret["response"]["debug"] = debug.str();
    ret["response"]["data"]["follow_sol_not_die"] = follow_sol_not_die;
    Json::Value v(Json::arrayValue);
    for (deque<int>::iterator iter = sol_not_die.begin(); iter != sol_not_die.end(); ++iter)
        v.append(*iter);
    ret["response"]["data"]["sol_not_die"] = v;
    Json::FastWriter writer;
    cout << writer.write(ret) << endl;
    exit(0);
}

bool dfs_die_in_step(const int step, const int round, Snake& snake) {
    if (!step) {
        ++cnt_sol_not_die;
        return false;
    }
    const Point ori_tail = snake.body.back();
    const bool grow = will_grow(round);
    if (!grow) snake.delete_back();
    bool any = false;
    for (int i = 0; i < 4; ++i) {
        if (!snake.valid_direction(i)) continue;
        const Point p = snake.body.front() + XY[i];
        snake.insert_front(p);
        if (!cnt_sol_not_die) sol_not_die[DIE_STEP-step] = i;
        if (cnt_sol_not_die < 2) {
            const bool die = dfs_die_in_step(step-1, round+1, snake);
            any = any || !die;
        }
        snake.delete_front();
    }
    if (!grow) snake.insert_back(ori_tail);
    return !any;
}

bool will_die_in_step(const int step, const int round, const Snake& ori) {
    memcpy(fmap, map, sizeof(map));
    sol_not_die.resize(DIE_STEP);
    Snake snake(ori.fork(fmap));
    bool res = dfs_die_in_step(step, round, snake);
    debug << "will die in " << step << " steps? " << (res ? "YES!" : "NO!") << endl;
    return res;
}

void board_bfs(const Point& source) {
    memset(board_dist, 0x3F, sizeof(board_dist));
    for (board_dist[source.x][source.y] = 0, Q.push(source); !Q.empty(); Q.pop()) {
        const Point &p = Q.front();
        const int d = board_dist[p.x][p.y];
        for (int i = 0; i < 4; ++i) {
            const Point t = p + XY[i];
            if (map[t.x][t.y] != '.' || board_dist[t.x][t.y] <= d+1) continue;
            board_dist[t.x][t.y] = d+1;
            Q.push(t);
        }
    }
}

int longest_path_direction(const Point& source) {
    memset(board_dist, 0x3F, sizeof(board_dist));
    int furthest_dist = 0;
    Point k(0, 0);
    for (board_dist[source.x][source.y] = 0, Q.push(source); !Q.empty(); Q.pop()) {
        const Point &p = Q.front();
        const int d = board_dist[p.x][p.y];
        if (d > furthest_dist) {
            furthest_dist = d;
            k = p;
        }
        for (int i = 0; i < 4; ++i) {
            const Point t = p + XY[i];
            if (map[t.x][t.y] != '.' || board_dist[t.x][t.y] <= d+1) continue;
            board_dist[t.x][t.y] = d+1;
            board_from[t.x][t.y] = i;
            Q.push(t);
        }
    }
    if (furthest_dist == 0) return -1;
    while (true) {
        const int d = board_from[k.x][k.y];
        const Point t = k - XY[d];
        if (t.x == source.x && t.y == source.y) return d;
        k = t;
    }
}

bool will_certainly_die_close(const int d, const Snake& ori) {
    static bool vis[MAXN][MAXN];
    memset(vis, 0, sizeof(vis));
    int cnt_empty = 0, earliest_leave = INF;
    const Point source = ori.body.front() + XY[d];
    for (Q.push(source), vis[source.x][source.y] = true; !Q.empty(); Q.pop()) {
        const Point u = Q.front();
        ++cnt_empty;
        for (int i = 0; i < 4; ++i) {
            const Point p = u + XY[i];
            if (map[p.x][p.y] == '0'+ori.id) {
                const int r = ori.body_leave_round(p);
                earliest_leave = std::min(earliest_leave, r);
            }
            if (map[p.x][p.y] != '.' || vis[p.x][p.y]) continue;
            vis[p.x][p.y] = true, Q.push(p);
        }
    }
    return earliest_leave-cur_round <= cnt_empty;
}

void decide() {
    const Point head(snake[0].body.front());
    const Point tail(snake[0].body.back());
    const bool grow = will_grow(cur_round);
    debug << (grow ? "will grow" : "won't grow") << endl;
    
    if (follow_sol_not_die) {
        debug << "follow sol_not_die" << endl;
        const int d = sol_not_die.front();
        sol_not_die.pop_front();
        follow_sol_not_die = sol_not_die.size();
        write_result(d);
    }

    debug << "trying longest_path_direction..." << endl;
    const int d = longest_path_direction(head);
    bool die = will_certainly_die_close(d, snake[0]);
    if (d != -1) {
        snake[0].insert_front(head+XY[d]);
        die = will_die_in_step(DIE_STEP, cur_round+1, snake[0]);
        snake[0].delete_front();
    }
    if (!die) {
        debug << "use longest_path_direction" << endl;
        follow_sol_not_die = cnt_sol_not_die == 1;
        write_result(d);
    }
    
    
    
    board_bfs(tail);
    int mindist = INF-1;
    vector<int> possible;
    debug << "trying board_dfs..." << endl;
    for (int i = 0; i < 4; ++i)
        if (snake[0].valid_direction(i)) {
            const Point p = head + XY[i];
            const int d = board_dist[p.x][p.y];
            if (d == 0 || d > mindist) continue;
            if (will_certainly_die_close(d, snake[0])) continue;
            
            snake[0].insert_front(p);
            debug << "testing direction " << i << ":  ";
            const bool die = will_die_in_step(DIE_STEP, cur_round+1, snake[0]);
            snake[0].delete_front();
            if (die) continue;
            
            if (d < mindist) {
                mindist = d;
                possible.clear();
            }
            if (d == mindist) possible.push_back(i);
        }
    if (!possible.empty()) {
        debug << "use board_bfs" << endl;
        const int d = rand()%possible.size();
        snake[0].insert_front(head+XY[d]);
        dfs_die_in_step(DIE_STEP, cur_round, snake[0]);
        follow_sol_not_die = cnt_sol_not_die == 1;
        write_result(possible[d]);
    }
    

    debug << "fallback to any valid direction" << endl;
    possible.clear();
    for (int i = 0; i < 4; ++i)
        if (snake[0].valid_direction(i))
            possible.push_back(i);
    write_result(possible.empty() ? -1 : possible[rand()%possible.size()]);
}

int main() {
    //freopen("input.json", "r", stdin);
    init();
    decide();
}
