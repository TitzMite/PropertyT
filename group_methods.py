
import datetime

def sorted_inverses(words):
    return sorted([w.inverse() for w in words])

def compute_products(left, right):
    prods = set()
    for l in left:
        for r in right:
            prods.add(l*r)
    return sorted(prods)

def cayley_ball(generators, radius):
    generators = list(generators)
    generators += [g.inverse() for g in generators]
    identity = generators[0].identity()
    ball = {identity}
    frontier = {identity}
    for _ in range(radius):
        new_frontier = set()
        for word in frontier:
            for generator in generators:
                new_word = word * generator
                if new_word not in ball:
                    ball.add(new_word)
                    new_frontier.add(new_word)
        frontier = new_frontier
    return sorted(ball)

def compute_table(left, right):
    print(datetime.datetime.now())
    print("starting to compute table")
    print(f"we have {len(left)} left words and {len(right)} right words")
    product_rows = []
    prods = set()
    total = len(left)
    for i, l in enumerate(left, start=1):
        row = []
        for r in right:
            product = l * r
            row.append(product)
            prods.add(product)
        product_rows.append(row)
        percent = 100 * i / total
        print(
            f"\rprogress computing table: {percent:.2f}%",
            end="",
            flush=True,
        )
    print()
    products = sorted(prods)
    print(f"we have {len(products)} products")
    index = {w: i for i, w in enumerate(products)}
    table = [
        [index[product] for product in row]
        for row in product_rows
    ]
    print("table created\n")
    return table