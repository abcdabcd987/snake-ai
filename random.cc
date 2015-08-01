#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <list>
#include <vector>
#include <iostream>

#include "jsoncpp/json.h"

const int MAXN = 25;
const int XY[4][2] = {{-1, 0}, {0, 1}, {1, 0}, {0, -1}};

struct Point {
    int x, y;
    Point(int a, int b): x(a), y(b) {}
};

int n, m;
char map[MAXN][MAXN+1];

struct Snake {
    int id;
    std::list<Point> body;
    bool inBody(int x, int y);
    void insertFront(const Point& p);
    void deleteBack(); 
    void move(int dir, int round); 
    void print() const; 
    bool validDirection(int dir);
} snake[2];

inline int getSnake(int x, int y) {
    if (map[x][y] == '0' || map[x][y] == '1') return map[x][y] - '0';
    return -1;
}

inline bool willGrow(int round) {
    if (round <= 9) return true;
    return (round-9) % 3 == 0;
}

inline bool Snake::inBody(int x, int y) {
    for (std::list<Point>::iterator iter = body.begin(); iter != body.end(); ++iter)
        if (x == iter->x && y == iter->y)
            return true;
    return false;
}

inline void Snake::insertFront(const Point& p) {
    body.push_front(p);
    map[p.x][p.y] = '0' + id;
}

inline void Snake::deleteBack() {
    const Point p(body.back());
    map[p.x][p.y] = '.';
    body.pop_back();
}

inline void Snake::move(int dir, int round) {
    const Point &p(body.front());
    insertFront(Point(p.x+XY[dir][0], p.y+XY[dir][1]));
    if (!willGrow(round)) deleteBack();
}

inline void Snake::print() const {
    for (std::list<Point>::const_iterator iter = body.begin(); iter != body.end(); ++iter)
        printf("(%d, %d) ", iter->x, iter->y);
    puts("");
}

inline bool Snake::validDirection(int dir) {
    const int x = body.front().x + XY[dir][0];
    const int y = body.front().y + XY[dir][1];
    if (x > n || y > m || x < 1 || y < 1) return false;
    if (map[x][y] != '.') return false;
    return true;
}

void printMap() {
    for (int i = 1; i <= n; ++i)
        printf("%s\n", map[i]);
}

void init() {
    memset(map, '.', sizeof(map));
    std::string str, tmp;
    while (getline(std::cin, tmp)) str += tmp;

    Json::Reader reader;
    Json::Value input;
    reader.parse(str, input);
    n = input["requests"][static_cast<Json::Value::UInt>(0)]["height"].asInt();
    m = input["requests"][static_cast<Json::Value::UInt>(0)]["width"].asInt();

    const int x = input["requests"][static_cast<Json::Value::UInt>(0)]["x"].asInt(); 
    snake[x].insertFront(Point(n, m));
    snake[1-x].insertFront(Point(1, 1));
    snake[0].id = 0, snake[1].id = 1;

    const int cnt_obstacle = input["requests"][static_cast<Json::Value::UInt>(0)]["obstacle"].size();
    for (int i = 0; i < cnt_obstacle; ++i) {
        const int x = input["requests"][static_cast<Json::Value::UInt>(0)]["obstacle"][static_cast<Json::Value::UInt>(0)]["x"].asInt();
        const int y = input["requests"][static_cast<Json::Value::UInt>(0)]["obstacle"][static_cast<Json::Value::UInt>(0)]["x"].asInt();
        map[x][y] = '#';
    }

    const int cnt_round = input["responses"].size();
    for (int i = 0; i < cnt_round; ++i) {
        const int dir0 = input["responses"][i]["direction"].asInt();
        const int dir1 = input["responses"][i+1]["direction"].asInt();
        snake[0].move(dir0, i);
        snake[1].move(dir1, i);
    }

    if (!willGrow(cnt_round)) {
        snake[0].deleteBack();
        snake[1].deleteBack();
    }

    srand(static_cast<unsigned>(time(0)) + cnt_round);
}

void decide() {
    std::vector<int> possible;
    for (int i = 0; i < 4; ++i)
        if (snake[0].validDirection(i))
            possible.push_back(i);

    Json::Value ret;
    ret["response"]["direction"] = possible.size() ? possible[rand()%possible.size()] : 0;
    Json::FastWriter writer;
    std::cout << writer.write(ret) << std::endl;
}

int main() {
    init();
    decide();
}
