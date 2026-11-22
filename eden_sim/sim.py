from typing import Literal, Set
import numpy as np
import pandas as pd
from pathlib import Path
from itertools import accumulate
from functools import reduce
import ast
from itertools import permutations

rng = np.random.default_rng(42)
CellState = Literal["edge", "center"]
CellAction = Literal["death", "birth", "skip"]


class Cell:
    def __init__(self, cell_id: int, parent_id: int | None,
                 state: CellState,
                 x: int, y: int, birth_date: int, generation: int,
                 death_rate: float, birth_rate: float,
                 mutations: Set[int]):
        self.cell_id = cell_id
        self.parent_id = parent_id
        self.state = state
        self.x, self.y = x, y
        self.birth_date = birth_date
        self.generation = generation
        self.death_rate = death_rate
        self.birth_rate = birth_rate
        self.mutations = mutations

    def choose_action(self) -> CellAction:
        """Choose an action based on cells birth, death probability. returns "birth", "death", or "skip". returning "birth"
        does not mean the cell multiplies, it may be in a state (surrounded state) where it can't reproduce (under certain models)"""
        if rng.random() < self.death_rate:
            return "death"
        elif rng.random() < self.birth_rate:
            return "birth"
        return "skip"


class World:
    def __init__(self, size: int,
                 time_step_hours: float = 2.0,
                 target_cells: int = 1000,
                 growth_model: Literal["boundary", "unrestricted"] = "boundary",
                 mutation_rate: float = 1.0,
                 log_file: str = None):

        # check if inputs make sense
        if size < 1 or target_cells < 1:
            raise ValueError("size and target_cells must be positive")

        if mutation_rate < 0 or mutation_rate > 1:
            raise ValueError(
                f"`mutation_rate` must be a float between [0, 1], but {mutation_rate} was given")

        if growth_model not in ("boundary", "unrestricted"):
            raise ValueError("growth model must be either boundary or unrestricted")

            # if the matrix is not enuf to even reach target_cells
        if target_cells > size * size:
            raise ValueError(
                "the matrix has too few positions compared to the upper limit set by target_cells")

        # if Path(log_file).is_file():
        #     raise ValueError("Log File Already Exists, Aborting")

        self.cell_count = 0
        self.size = size
        self.grid = np.full((size, size), None, dtype = Cell)

        # each cell is given a unique identifying ID, that starts from 1 and goes up
        self.next_id = 1

        # to track sim time --> can print out for extra detail if needed/want
        self.step = 0
        self.time_step_hours = time_step_hours

        # to know when to stop sim (default being 1000 like in the paper)
        self.target_cells = target_cells

        # type of model to use --> boundary means don't push to make space when the 8 surrounding cells are filled, unrestricted means do that
        self.growth_model = growth_model

        # mutations follows a Jukes-Cantor Model (equal nucleotide frequency, independent mutation rate at each position, equal substitution rates between any 2 nucleotides)
        # and an infinite site assumption (a specific base is only ever mutated once -> one strain cannot mutate a base twice AND 2 strains cannot independently mutate the same base
        # to implement infinite site assumption, whenever we create a mutation, we increase the size of the genome by 1, and mutate it (so it's impossible to mutate a previously added base)
        # mutations are stored in a set of bases i.e. {1, 5, 8} = this individual has a mutation in position 1, 5, 8
        #
        # next base position to add
        self.next_mutation_position = 1
        # probability of offspring to get a mutation (1.0 in the paper)
        self.mutation_rate = mutation_rate

        # Path where logs will be generated as CSV
        self.log_file = Path(log_file)
        with open(self.log_file, "w") as file:
            file.write("id\tparent_id\tbirth_date\tdeath_date\tmutations\n")

    def update(self):
        # after adding remove cell function (which happens when cell action is death), sometimes sim keeps running when there's already nothing
        # repetitive bc when there's nothing, the next step is also going to be nothing
        # stop returning if the above case or if it goes past the set upper limit on # of cells
        if self.cell_count >= self.target_cells or self.cell_count == 0:
            return

        # update self.step and self.time_hours accordingly
        self.step += 1

        # Every cell that exists at the start of update should be updated, in random order to not be biased
        cells = self._get_all_cells_shuffled()

        for cell in cells:
            action = cell.choose_action()

            match action:
                case "death":
                    self.remove_cell(cell)
                case "birth":
                    self.birth_new_cell(cell)
            # print(f"cell {cell.cell_id} chooses to {action}")

            if self.cell_count >= self.target_cells:
                break

    def _get_all_cells_shuffled(self):
        """Return a shuffled list of all currently alive cells"""
        cells = self.grid.flatten()
        cells = cells[cells != None]
        rng.shuffle(cells)
        return cells

    def _set_cell_at_position(self, cell: Cell):
        self.grid[cell.x, cell.y] = cell
        self.cell_count += 1

    def _get_next_cell_id(self):
        """return next unused id number and increment counter"""
        self.next_id += 1
        return self.next_id - 1

    def birth_new_cell(self, parent: Cell):
        """Look for an empty neighboring position from the parent cell and create a new Cell there"""
        empty_neighboring_indices = []
        for x in range(parent.x - 1, parent.x + 2):
            for y in range(parent.y - 1, parent.y + 2):
                if x == parent.x and y == parent.y:
                    continue
                if x < 0 or y < 0 or x >= self.size or y >= self.size:
                    continue
                if self.grid[x, y] is None:
                    empty_neighboring_indices.append((x, y))

        if empty_neighboring_indices:
            cell_position = empty_neighboring_indices[rng.integers(len(empty_neighboring_indices))]
        elif self.growth_model == "boundary":
            # if boundary, u can't do anything abt it
            return
        else:
            # if unrestricted, make space by pushing surroundings
            cell_position = self.push_to_make_space(parent)

            if cell_position is None:
                return

        # mutations are inherited from the parent, and if cell mutates in process, add a new mutation
        mutations = parent.mutations.copy()
        if rng.random() < self.mutation_rate:
            mutations.add(self.next_mutation_position)
            self.next_mutation_position += 1

        self._set_cell_at_position(Cell(
            cell_id = self._get_next_cell_id(), parent_id = parent.cell_id,
            state = "edge",
            x = cell_position[0], y = cell_position[1],
            birth_date = self.step,
            generation = parent.generation + 1,
            death_rate = parent.death_rate, birth_rate = parent.birth_rate,
            mutations = mutations
        ))

    def seed_cell(self, x: int, y: int, birth_rate: float, death_rate: float):
        """Used at the start of the simulation to seed the founding cell"""
        self._set_cell_at_position(
            Cell(cell_id = self._get_next_cell_id(), parent_id = None, state = "edge",
                 x = x, y = y, birth_date = self.step, death_rate = death_rate,
                 generation = 1,
                 birth_rate = birth_rate,
                 mutations = set())
        )

    def remove_cell(self, cell: Cell):
        # On death, add a log of that cell (this will need to be changed if we want to store how a cell's position changes through time)
        with open(self.log_file, "a") as file:
            file.write(
                f"{cell.cell_id}\t{cell.parent_id}\t{cell.birth_date}\t{self.step}\t{cell.mutations}\n")

        self.grid[cell.x, cell.y] = None
        self.cell_count -= 1

    def push_to_make_space(self, parent: Cell) -> tuple[int, int] | None:
        # output of None meaning it failed to make space --> reached matrix limits
        # 8 lattice directions, including diagonals
        directions = [
            (dx, dy)
            for dx in (-1, 0, 1)
            for dy in (-1, 0, 1)
            if (dx, dy) != (0, 0)
        ]

        dx, dy = directions[rng.integers(len(directions))]

        x = parent.x + dx
        y = parent.y + dy

        occupied_positions = []

        # find first empty position in chosen direction
        while True:
            if x < 0 or y < 0 or x >= self.size or y >= self.size:
                # raise RuntimeError("pushing reached grid boundary, u can try again w/ larger grid inputs")
                return None

            if self.grid[x, y] is None:
                break

            occupied_positions.append((x, y))

            x += dx
            y += dy

        # do in reverse to move farthest cell first so that everything moves cleanly/no overlap
        for old_x, old_y in reversed(occupied_positions):
            cell = self.grid[old_x, old_y]

            new_x = old_x + dx
            new_y = old_y + dy

            self.grid[new_x, new_y] = cell
            self.grid[old_x, old_y] = None

            cell.x = new_x
            cell.y = new_y

        # this neighboring position now empty and birth_new_cell() will put the cell in this empty slot
        return parent.x + dx, parent.y + dy

    def generate_strands(self):
        """
        Uses information in log file and appends a generated DNA sequence to it. Updates the log file
        :return:
        Dataframe of log with new column.
        """
        df = pd.read_csv(
            self.log_file, sep = '\t',
            dtype = {"id": "string", "parent_id": "string", "birth_date": "int",
                     "death_date": "int"}
        )
        df['mutations'] = df['mutations'].apply(ast.literal_eval)
        seq_length = max(
            (max(x, default = 0) for x in df["mutations"]),
            default = 0,
        )

        bases = ["A", "T", "G", "C"]
        perms = list(permutations(bases, r = 2))
        mapping = {
            i: perms[rng.integers(0, len(perms))]
            for i in range(seq_length)
        }

        df['seq'] = df['mutations'].map(
            lambda mut:
            "".join(mapping[i][1] if i in mut else mapping[i][0] for i in range(seq_length))
        )

        # Write result to logfile
        df.to_csv(self.log_file, sep = "\t", index = False)
        return df

    def end_sim(self):
        """Kill all remaining alive cells and log them"""
        cells = self._get_all_cells_shuffled()
        for cell in cells:
            self.remove_cell(cell)

        self.generate_strands()


# Helper Functions
def visualize_world(world: World):
    """hacky way to visualize each world grid"""
    grid = world.grid
    for row in grid:
        print("  ".join(["..." if x is None else f"{x.cell_id:03d}" for x in row]), flush = True)
