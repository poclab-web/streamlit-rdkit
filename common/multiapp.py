"""Frameworks for running multiple Streamlit applications as a single app.
"""
import streamlit as st

class MultiApp:
    """Framework for combining multiple streamlit applications.
    Usage:
        def foo():
            st.title("Hello Foo")
        def bar():
            st.title("Hello Bar")
        app = MultiApp()
        app.add_app("Foo", foo)
        app.add_app("Bar", bar)
        app.run()
    It is also possible keep each application in a separate file.
        import foo
        import bar
        app = MultiApp()
        app.add_app("Foo", foo.app)
        app.add_app("Bar", bar.app)
        app.run()
    """
    def __init__(self):
        self.apps = []

    def add_app(self, title, func, category=None, show_code=True):
        """Adds a new application.
        Parameters
        ----------
        title : str
            Title of the app. Appears in the dropdown in the sidebar.
        func : callable
            The Python function to render this app.
        category : str, optional
            The category under which the app appears.
        show_code : bool, optional
            Whether to show the source code for this app. Defaults to True.
        """
        self.apps.append({
            "title": title,
            "function": func,
            "category": category,
            "show_code": show_code  # 新しいオプション
        })

    def run(self):
        categorized_apps = [app for app in self.apps if app['category'] is not None]
        uncategorized_apps = [app for app in self.apps if app['category'] is None]

        with st.sidebar:
            st.title("ページ選択")

            # カテゴリーなしのアプリを表示
            if uncategorized_apps:
                st.markdown("### その他")
                for app in uncategorized_apps:
                    if st.button(app['title'], key=app['title']):
                        self.display_app(app)  # 修正: アプリ表示を専用メソッドに分離
                        return

            # カテゴリーごとの整理
            categories = sorted(set(app['category'] for app in categorized_apps))
            for category in categories:
                with st.expander(category):  # 各カテゴリをエクスパンダーで整理
                    selected_title = st.radio(
                        "選択してください:",
                        [app['title'] for app in categorized_apps if app['category'] == category],
                        key=category
                    )
                    for app in categorized_apps:
                        if app['title'] == selected_title:
                            self.display_app(app)  # 修正: アプリ表示を専用メソッドに分離
                            return

    def display_app(self, app):
        """Displays the selected app and optionally its source code.
        Parameters
        ----------
        app : dict
            The app to display, including its function and metadata.
        """
        st.title(f"現在のアプリ: {app['title']}")


        # アプリを実行
        try:
            app['function']()
        except Exception as e:
            st.error(f"アプリの実行中にエラーが発生しました: {e}")

