# 化学方程式精确配平 API

基于 FastAPI 的化学方程式配平服务。提交反应物与生成物分子式后，服务端构建
**含电荷行的元素守恒矩阵**，用**有理数高斯消元**（`fractions.Fraction`，全程精确、
无浮点误差）求零空间：

- **零空间恰为一维**且基向量可化为**全正系数** → 返回按最大公约数约简的
  **最小正整数系数**，并附逐元素、逐电荷的守恒核算；
- 否则返回 **HTTP 422**，说明零空间维数（`no_solution` / `underdetermined`）
  或无法全正的原因（`not_all_positive`）。

## 快速启动（单条命令）

```bash
docker compose up --build
```

依赖（`requirements.txt`）在容器**启动时**由 `start.sh` 自动安装，随后拉起 uvicorn。

启动后访问：

| 地址 | 说明 |
| --- | --- |
| http://localhost:8000/docs | Swagger 交互式文档 |
| http://localhost:8000/redoc | ReDoc 文档 |
| http://localhost:8000/health | 健康检查 |

## 配置

通过环境变量配置（可在 `docker-compose.yml` 的 `environment` 中修改）：

| 变量 | 默认值 | 说明 |
| --- | --- | --- |
| `PORT` | `8000` | 容器内 uvicorn 监听端口（Compose 已映射到宿主机 8000） |
| `HOST` | `0.0.0.0` | 监听地址 |

## 分子式语法

- 元素与下标：`H2O`、`Fe2O3`
- 嵌套括号（`()` `[]` `{}`，可任意嵌套）：`Ca3(PO4)2`、`K4[Fe(CN)6]`
- 结晶水（`.`、`·`、`•` 或 `*` 分隔，段首可带系数）：`CuSO4·5H2O`、`Na2CO3.10H2O`
- 末尾离子电荷（`^` + 数字 + 正负号）：`SO4^2-`、`Fe^3+`、`H^+`、`OH^-`、`Al^{3+}`

## API

### `POST /api/v1/balance` — 配平

```bash
curl -s -X POST http://localhost:8000/api/v1/balance \
  -H 'Content-Type: application/json' \
  -d '{"reactants": ["KMnO4", "HCl"], "products": ["KCl", "MnCl2", "H2O", "Cl2"]}'
```

响应（200）：

```json
{
  "equation": "2 KMnO4 + 16 HCl -> 2 KCl + 2 MnCl2 + 8 H2O + 5 Cl2",
  "coefficients": {"reactants": [2, 16], "products": [2, 2, 8, 5]},
  "null_space_dimension": 1,
  "verification": {
    "elements": [
      {"element": "K", "reactants": 2, "products": 2, "balanced": true}
    ],
    "charge": {"reactants": 0, "products": 0, "balanced": true},
    "balanced": true
  }
}
```

离子方程式示例（电荷行参与守恒）：

```bash
curl -s -X POST http://localhost:8000/api/v1/balance \
  -H 'Content-Type: application/json' \
  -d '{"reactants": ["MnO4^-", "H^+", "Fe^2+"], "products": ["Mn^2+", "Fe^3+", "H2O"]}'
# MnO4^- + 8 H^+ + 5 Fe^2+ -> Mn^2+ + 5 Fe^3+ + 4 H2O   （两侧电荷均为 +17）
```

无法唯一配平时返回 422：

```json
{
  "detail": {
    "error": "underdetermined",
    "message": "The null space of the conservation matrix has dimension 2: ...",
    "null_space_dimension": 2
  }
}
```

`error` 取值：

| error | 含义 |
| --- | --- |
| `no_solution` | 零空间维数为 0，只有全零平凡解，无法配平 |
| `underdetermined` | 零空间维数 > 1，存在无穷多组解，无唯一系数 |
| `not_all_positive` | 零空间一维，但含零或异号系数，无法化为全正 |
| `invalid_formula` | 分子式语法错误（未知元素、括号不配对、电荷格式错误等） |

### `POST /api/v1/validate` — 系数校验

校验一组给定系数是否满足元素与电荷守恒（始终返回 200，看 `balanced` 字段）：

```bash
curl -s -X POST http://localhost:8000/api/v1/validate \
  -H 'Content-Type: application/json' \
  -d '{"reactants": ["H2", "O2"], "products": ["H2O"],
       "coefficients": {"reactants": [2, 1], "products": [2]}}'
```

## 本地开发与测试

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements-dev.txt
pytest tests/                 # 全部测试（核心算法 + API）
python -m unittest tests.test_chemistry   # 仅核心算法（纯标准库，无需安装依赖）
uvicorn app.main:app --reload
```

## 项目结构

```
app/
  chemistry.py   # 分子式解析、守恒矩阵、有理数高斯消元求零空间、配平与核算
  schemas.py     # Pydantic 请求/响应模型
  main.py        # FastAPI 应用与路由
tests/           # 核心算法与 API 测试
Dockerfile       # 镜像（依赖随容器启动安装，见 start.sh）
docker-compose.yml
start.sh         # 启动脚本：pip install -> uvicorn
```
