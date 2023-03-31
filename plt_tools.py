# %% [code]

def get_figsize(
    columnwidth=4, wf=1.0, hf_rel=(5.0 ** 0.5 - 1.0) / 2.0, hf_abs=None, unit="inch"
):
    # function to be used with
    # guidelines from elsevier for sizing of artwork
    # https://www.elsevier.com/authors/policies-and-guidelines/artwork-and-media-instructions/artwork-sizing
    """A função aceita uma série de argumentos para a precisa customização da
    largura e altura de uma figura que será produzida com Matplotlib.
    
    Maiores detalhes em:
    www.fschuch/blog/2020/10/14/graficos-com-qualidade-de-publicacao-em-python-com-matplotlib
    
    © 2020 Felipe N. Schuch, sob os termos da CC BY SA 4.0.
    (www.creativecommons.org/licenses/by-sa/4.0)
    
    Parâmetros
    ----------
    columnwidth : float
        Largura da página ou da coluna de texto (o padrão é 4)
        Obtenha isso em LaTeX usando \the\columnwidth
    wf : float
        A fração da largura que será utilizada pela figura (o padrão é 1.0)
    hf_rel : float
        Altura da figura, valor relativo em relação à largura (o padrão é a proporção áurea)
    hf_abs : float
        Altura da figura em termos absolutos, caso desejado (o padrão é None)
        Obtenha isso em LaTeX usando \the\textheight
    unit : float
        Unidade de comprimento utilizada para `columnwidth` e `hf_abs`,
        as opções suportadas são "inch" (polegada), "mm", "cm" e "pt" (Pontos tipográfico)
        (o padrão é "inch")
    Retorno
    -------
    set
        Retorna um set contendo a largura e a altura especificada para a figura.
    Apura
    -------
    ValueError
        Caso a unidade não seja suportada pela função.
    Exemplos
    -------
    
    >>> import matplotlib.pyplot as plt
    
    >>> get_figsize()
    (4.0, 2.4721359549995796)
    >>> get_figsize(columnwidth=16, unit='cm', hf_abs=9)
    (6.299212598425196, 3.543307086614173)
    >>> get_figsize(columnwidth=4, unit='inch', hf_rel = 1.0)
    (4.0, 4.0)

    >>> plt.rcParams.update({
    ... 'figure.figsize' : get_figsize(columnwidth=160, wf=0.75, unit='mm', hf_rel=1)
    ... })
    """

    # Dessa maneira, unit não será sensível a letras maiúsculas e minúsculas
    unit = unit.lower()

    # Converte unidades para polegadas, conforme esperado por Matplotlib
    conversion = dict(inch=1.0, mm=25.4, cm=2.54, pt=72.0,)

    if unit in conversion.keys():
        fig_width = columnwidth / conversion[unit]
        if hf_abs is not None:
            fig_height = hf_abs / conversion[unit]
    else:
        raise ValueError(f"unit deve ser: {conversion.keys()}")

    # A figura será apenas uma fração da largura útil da página
    fig_width *= wf

    # Caso hf_abs não seja definido, a altura será uma fração da largura
    if hf_abs is None:
        fig_height = fig_width * hf_rel

    # Retorna a largura e altura especificada para a figura
    return (fig_width, fig_height)

def set_plt_rc(plt, for_paper = False):
    if for_paper:
        ## https://www.elsevier.com/authors/policies-and-guidelines/artwork-and-media-instructions/artwork-sizing
        plt.rc('font', family='serif', size=7)#, serif= ['Computer Modern']) # controls default text sizes
    else:
        plt.rc('font', family='serif', size=12) # controls default text sizes
        plt.rc('axes', titlesize=14)     # fontsize of the axes title
        plt.rc('axes', labelsize=13)    # fontsize of the x and y labels
        plt.rc('xtick', labelsize=13)    # fontsize of the tick labels
        plt.rc('ytick', labelsize=13)    # fontsize of the tick labels
        plt.rc('legend', fontsize=12)    # legend fontsize
        plt.rc('figure', titlesize=14)  # fontsize of the figure title
        
    return plt