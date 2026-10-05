"""Bundled, searchable documentation; available before scientific data load."""
from importlib.resources import files
import re
from . import __version__

DOCUMENTS={'User manual':'USER_MANUAL.md','Installation guide (README)':'README.md'}

def read_document(name):
    # Resource names are a fixed allowlist, never a path supplied by a visitor.
    return files('microc_explorer').joinpath('assets',DOCUMENTS[name]).read_text(encoding='utf-8')

def sections(markdown):
    parts=re.split(r'(?m)^## ',markdown)
    result=[]
    for part in parts[1:]:
        title,_,body=part.partition('\n')
        if title.strip()!='Contents':
            result.append((title.strip(),body.strip()))
    return result

def render_manual():
    import streamlit as st
    st.header('User manual')
    st.caption(f'Micro-C Explorer {__version__} · usage, installation and administration')
    name=st.selectbox('Documentation',list(DOCUMENTS),key='manual_document')
    text=read_document(name)
    query=st.text_input('Search documentation',placeholder='For example: gene, Docker, memory, 2 Mb',
                        max_chars=120,key='manual_search').strip().casefold()
    all_sections=sections(text)
    topic=st.selectbox('Topic',['All topics']+[t for t,_ in all_sections],key=f'manual_topic_{name}')
    selected=[(title,body) for title,body in all_sections
              if (topic=='All topics' or title==topic) and
              (not query or query in (title+'\n'+body).casefold())]
    if not selected:
        st.info('No matching topics. Clear the search or choose All topics.')
    else:
        st.caption(f'{len(selected)} topic(s) shown')
        for i,(title,body) in enumerate(selected):
            with st.expander(title,expanded=bool(query) or topic!='All topics' or i==0):
                st.markdown(body)
    left,right=st.columns(2)
    left.download_button('Download user manual',read_document('User manual'),
                         'MicroC_Explorer_User_Manual.md','text/markdown',key='manual_download')
    right.download_button('Download README',read_document('Installation guide (README)'),
                          'README.md','text/markdown',key='readme_download')
