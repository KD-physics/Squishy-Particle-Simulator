"""
Color palettes for particle rendering. Apply via System.set_color_palette(name)
or directly via apply_palette(sys, name).

Palettes are lists of hex color strings (matplotlib-compatible). Each palette
has a thematic group; pick by name from PALETTES.
"""
import random as _random


PALETTES = {
    # 70's pastel — dull apple red, peach, butter, mint, purple, cyan
    'palette1':  ['#cc5555', '#ffb347', '#ffff80', '#77dd77', '#b16ad9', '#00e1e1'],

    # Muted modern — teal, pale sand, slate blue, coral
    'palette2':  ['#008080', '#ece0c8', '#708090', '#f08080'],

    # Earthy — soft olive, warm beige, terracotta, cool gray
    'palette3':  ['#aabd91', '#ded2ba', '#c68667', '#b0b6ba'],

    # Cool muted — blue-gray, blush, taupe, navy
    'palette4':  ['#8faab4', '#f4dddc', '#c4b1a2', '#3c4c59'],

    # Bright warm — tomato, orange, yellow, pale green, lemon, lime
    'palette5':  ['#ff6347', '#ffa500', '#ffdf00', '#c0f1c0', '#fffacd', '#32cd32'],

    # Bright cool — sky blue, dodger, hot pink, lavender, azure, red-orange
    'palette6':  ['#00bfff', '#1e7cf5', '#ff69b4', '#ebd7e6', '#dcebeb', '#ff4500'],

    # Vibrant magenta/green — magenta, spring green, lime, lavender, indigo
    'palette7':  ['#ff6eff', '#00ff7f', '#adff2f', '#ffebfa', '#9b50d2'],

    # Bold neon — deep pink, cyan, violet, lawn green
    'palette8':  ['#ff1493', '#00ffff', '#ee82ee', '#7cfc00'],

    # Winter — powder blue, light blue, lemon, sandy, azure, steel blue
    'palette9':  ['#b0e0e6', '#add8e6', '#fffacd', '#ea9a56', '#f0ffff', '#4682b4'],

    # Spring — light pink, pale green ×2, light blue, lavender, deep blue
    'palette10': ['#ffa2ad', '#8efb8e', '#d4fbd4', '#add8e6', '#ffebfa', '#5a96c8'],

    # Summer — orange, tomato, olive, lemon, sky
    'palette11': ['#ffa500', '#ff6347', '#bab86c', '#fffacd', '#87ceeb'],

    # Fall — indian red, gold, coffee, chocolate, sandy, moccasin
    'palette12': ['#cd5c5c', '#d5b56e', '#9c6f44', '#d2691e', '#ea9a56', '#f5daab'],
}


def apply_palette(sys_h, name, seed=None):
    """Randomly assign one color per particle from the named palette.

    Writes into sys_h._particle_colors. Particles already given an explicit
    color via set_driven_particles or _particle_colors[i]= are overwritten,
    so call this BEFORE coloring driven/pinned rings.
    """
    if name not in PALETTES:
        raise ValueError(f"unknown palette {name!r}; available: {sorted(PALETTES)}")
    colors = PALETTES[name]
    rng = _random.Random(seed)
    n = len(sys_h._particles)
    sys_h._particle_colors = [rng.choice(colors) for _ in range(n)]
    return sys_h._particle_colors
